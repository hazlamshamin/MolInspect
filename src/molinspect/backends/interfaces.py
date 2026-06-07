"""Optional biological-interface backends.

PISA is treated as authoritative biological-interface evidence when available.
MolInspect's cheaper inter-chain contact interface remains a separate heuristic
landmark so callers can tell contact geometry from PISA surface/assembly output.
"""

from __future__ import annotations

import os
import re
import shutil
import subprocess
import tempfile
import xml.etree.ElementTree as ET
from dataclasses import dataclass
from functools import lru_cache
from pathlib import Path
from typing import Any

from .common import short_process_error, source_cache_key, static_structure_source
from ..objects import object_id_for_residue
from ..relations import chain_for_residue
from ..schemas import ObjectRef
from ..definitions import definition

PISA_BACKEND = "PISA"
PISA_TIMEOUT_S = 120


@dataclass(frozen=True, slots=True)
class BiologicalInterfaceRecord:
    """PISA interface result mapped to residue selections."""

    id: str
    name: str
    backend: str
    pisa_id: str
    chains: tuple[str, ...]
    pisa_chain_ids: tuple[str, ...]
    interface_area_A2: float | None
    solvation_energy_kcal_mol: float | None
    stabilization_energy_kcal_mol: float | None
    pvalue: float | None
    residue_ids: tuple[str, ...]
    residue_expressions: tuple[str, ...]
    annotations: dict[str, Any]


@dataclass(frozen=True, slots=True)
class BiologicalInterfaceStore:
    """Detected biological interfaces and their residue memberships."""

    backend: str | None = None
    records: tuple[BiologicalInterfaceRecord, ...] = ()
    limitations: tuple[str, ...] = ()


def detect_biological_interfaces(
    universe: Any,
    source_files: tuple[Path, ...],
) -> BiologicalInterfaceStore:
    """Run PISA when available and map reported interfaces onto loaded residues."""

    source = static_structure_source(source_files)
    if source is None:
        return BiologicalInterfaceStore(
            limitations=("PISA biological-interface detection requires a static PDB/mmCIF source.",)
        )

    residues_by_key = _residues_by_chain_resid(universe)
    store = _cached_pisa_interfaces(source_cache_key(source))
    return _hydrate_biological_interface_store(store, residues_by_key)


def pisa_is_available() -> bool:
    """Return whether a local PISA command is discoverable."""

    return _pisa_command() is not None


def object_ref_for_biological_interface(record: BiologicalInterfaceRecord) -> ObjectRef:
    """Return a public object reference for one PISA biological interface."""

    return ObjectRef(
        id=record.id,
        type="biological_interface",
        name=record.name,
        annotations=record.annotations,
    )


def selection_expression_for_biological_interface(record: BiologicalInterfaceRecord) -> str:
    """Return an MDAnalysis selection expression for PISA interface residues."""

    return " or ".join(f"({expression})" for expression in record.residue_expressions)


@lru_cache(maxsize=16)
def _cached_pisa_interfaces(source_key: tuple[str, int, int]) -> BiologicalInterfaceStore:
    command = _pisa_command()
    if command is None:
        return BiologicalInterfaceStore(
            limitations=(
                "PISA is unavailable. Install CCP4/PISA and expose the `pisa` command, "
                "or set MOLINSPECT_PISA_COMMAND, to enable biological-interface objects.",
            )
        )

    source = Path(source_key[0])
    return _run_pisa(source, command)


def _run_pisa(source: Path, command: str) -> BiologicalInterfaceStore:
    with tempfile.TemporaryDirectory(prefix="molinspect-pisa-") as temp_dir:
        temp_path = Path(temp_dir)
        config_args, config_limitations = _pisa_config_args(temp_path)
        session_name = _pisa_session_name(source)
        analyse_cmd = [command, session_name, "-analyse", str(source), *config_args]
        try:
            completed = subprocess.run(
                analyse_cmd,
                cwd=temp_path,
                check=False,
                capture_output=True,
                text=True,
                timeout=PISA_TIMEOUT_S,
            )
        except subprocess.TimeoutExpired:
            return BiologicalInterfaceStore(
                backend=PISA_BACKEND,
                limitations=(f"PISA timed out during analysis after {PISA_TIMEOUT_S} seconds.",),
            )
        if completed.returncode != 0:
            detail = short_process_error(completed.stderr or completed.stdout)
            return BiologicalInterfaceStore(
                backend=PISA_BACKEND,
                limitations=tuple(
                    dict.fromkeys(
                        (
                            *config_limitations,
                            f"PISA analysis failed: {detail}",
                        )
                    )
                ),
            )

        xml_result = _pisa_interfaces_xml(command, session_name, config_args, temp_path)
        _erase_pisa_session(command, session_name, config_args, temp_path)
        if isinstance(xml_result, BiologicalInterfaceStore):
            return xml_result
        xml_path, provenance = xml_result
        try:
            records = _parse_pisa_interfaces_xml(xml_path)
        except Exception as exc:
            return BiologicalInterfaceStore(
                backend=PISA_BACKEND,
                limitations=(f"PISA interfaces XML could not be parsed: {exc}",),
            )
        if not records:
            return BiologicalInterfaceStore(
                backend=PISA_BACKEND,
                limitations=("PISA completed but reported no usable biological interfaces.",),
            )
        hydrated_records = tuple(
            _record_with_input_provenance(record, provenance) for record in records
        )
        return BiologicalInterfaceStore(backend=PISA_BACKEND, records=hydrated_records)


def _pisa_interfaces_xml(
    command: str,
    session_name: str,
    config_args: list[str],
    cwd: Path,
) -> tuple[Path, str] | BiologicalInterfaceStore:
    xml_path = cwd / "interfaces.xml"
    commands = (
        ([command, session_name, "-xml", "interfaces", *config_args], "PISA XML interfaces report"),
        ([command, session_name, "-xml", *config_args], "PISA XML report"),
    )
    errors: list[str] = []
    for xml_cmd, provenance in commands:
        try:
            completed = subprocess.run(
                xml_cmd,
                cwd=cwd,
                check=False,
                capture_output=True,
                text=True,
                timeout=PISA_TIMEOUT_S,
            )
        except subprocess.TimeoutExpired:
            errors.append(f"{xml_cmd[2]} timed out after {PISA_TIMEOUT_S} seconds")
            continue
        if completed.returncode != 0:
            errors.append(short_process_error(completed.stderr or completed.stdout))
            continue
        if not completed.stdout.strip():
            errors.append("PISA XML command produced empty stdout")
            continue
        xml_path.write_text(completed.stdout)
        return xml_path, provenance
    return BiologicalInterfaceStore(
        backend=PISA_BACKEND,
        limitations=(f"PISA XML retrieval failed: {' | '.join(errors)[:600]}",),
    )


def _erase_pisa_session(command: str, session_name: str, config_args: list[str], cwd: Path) -> None:
    try:
        subprocess.run(
            [command, session_name, "-erase", *config_args],
            cwd=cwd,
            check=False,
            capture_output=True,
            text=True,
            timeout=10,
        )
    except Exception:
        return


def _parse_pisa_interfaces_xml(path: Path) -> tuple[BiologicalInterfaceRecord, ...]:
    """Parse PISA interfaces XML into backend-neutral records."""

    root = ET.parse(path).getroot()
    if root.tag == "INTERFACE" or root.find(".//INTERFACE") is not None:
        return _parse_pisa_legacy_interface_xml(root)

    records: list[BiologicalInterfaceRecord] = []
    for index, element in enumerate(root.findall(".//interface"), start=1):
        pisa_id = _text(element, "id") or str(index)
        residue_entries = _interface_residue_entries(element)
        if not residue_entries:
            continue
        chains = tuple(dict.fromkeys(chain for chain, _, _ in residue_entries))
        pisa_chain_ids = tuple(
            dict.fromkeys(
                chain_id
                for molecule in element.findall("molecule")
                if (chain_id := _text(molecule, "chain_id")) is not None
            )
        )
        interface_definition = definition("biological_interface")
        interface_id = f"biological_interface:pisa:{_id_part(pisa_id)}"
        interface_area = _float_text(element, "int_area")
        solvation_energy = _float_text(element, "int_solv_en")
        stabilization_energy = _float_text(element, "stab_en")
        pvalue = _float_text(element, "pvalue")
        residue_ids = tuple(
            f"residue:{_id_part(chain)}:{_id_part(resid)}:{_id_part(resname)}"
            for chain, resid, resname in residue_entries
        )
        residue_expressions = tuple(
            _residue_expression_from_parts(chain, resid, resname)
            for chain, resid, resname in residue_entries
        )
        annotations: dict[str, Any] = {
            "definition_id": "biological_interface",
            "definition_source": interface_definition.source,
            "reference_keys": list(interface_definition.reference_keys),
            "backend": PISA_BACKEND,
            "method": interface_definition.parameters["method"],
            "pisa_interface_id": pisa_id,
            "chains": list(chains),
            "pisa_chain_ids": list(pisa_chain_ids),
            "interface_area_A2": interface_area,
            "solvation_energy_kcal_mol": solvation_energy,
            "stabilization_energy_kcal_mol": stabilization_energy,
            "pvalue": pvalue,
            "bond_counts": _bond_counts(element),
            "participating_residues": list(dict.fromkeys(residue_ids)),
            "participating_residue_count": len(set(residue_ids)),
            "limitation": interface_definition.limitations[0],
        }
        records.append(
            BiologicalInterfaceRecord(
                id=interface_id,
                name=f"pisa_biological_interface {pisa_id}",
                backend=PISA_BACKEND,
                pisa_id=pisa_id,
                chains=chains,
                pisa_chain_ids=pisa_chain_ids,
                interface_area_A2=interface_area,
                solvation_energy_kcal_mol=solvation_energy,
                stabilization_energy_kcal_mol=stabilization_energy,
                pvalue=pvalue,
                residue_ids=tuple(dict.fromkeys(residue_ids)),
                residue_expressions=tuple(dict.fromkeys(residue_expressions)),
                annotations=annotations,
            )
        )
    return tuple(records)


def _parse_pisa_legacy_interface_xml(root: ET.Element) -> tuple[BiologicalInterfaceRecord, ...]:
    """Parse legacy PDBePISA single-interface XML into records."""

    elements = [root] if root.tag == "INTERFACE" else root.findall(".//INTERFACE")
    records: list[BiologicalInterfaceRecord] = []
    for index, element in enumerate(elements, start=1):
        pisa_id = _text(element, "INTERFACENO") or str(index)
        residue_entries = _legacy_interface_residue_entries(element)
        if not residue_entries:
            continue
        chains = tuple(dict.fromkeys(chain for chain, _, _ in residue_entries))
        interface_definition = definition("biological_interface")
        interface_id = f"biological_interface:pisa:{_id_part(pisa_id)}"
        interface_area = _legacy_interface_area(element)
        solvation_energy = _legacy_solvation_energy(element)
        pvalue = _legacy_pvalue(element)
        residue_ids = tuple(
            f"residue:{_id_part(chain)}:{_id_part(resid)}:{_id_part(resname)}"
            for chain, resid, resname in residue_entries
        )
        residue_expressions = tuple(
            _residue_expression_from_parts(chain, resid, resname)
            for chain, resid, resname in residue_entries
        )
        annotations: dict[str, Any] = {
            "definition_id": "biological_interface",
            "definition_source": interface_definition.source,
            "reference_keys": list(interface_definition.reference_keys),
            "backend": PISA_BACKEND,
            "method": interface_definition.parameters["method"],
            "pisa_interface_id": pisa_id,
            "chains": list(chains),
            "pisa_chain_ids": list(chains),
            "interface_area_A2": interface_area,
            "interface_area_by_partner_A2": _legacy_partner_areas(element),
            "solvation_energy_kcal_mol": solvation_energy,
            "stabilization_energy_kcal_mol": None,
            "pvalue": pvalue,
            "bond_counts": _legacy_bond_counts(element),
            "participating_residues": list(dict.fromkeys(residue_ids)),
            "participating_residue_count": len(set(residue_ids)),
            "limitation": interface_definition.limitations[0],
            "xml_variant": "pdbe_pisa_legacy_interface",
        }
        records.append(
            BiologicalInterfaceRecord(
                id=interface_id,
                name=f"pisa_biological_interface {pisa_id}",
                backend=PISA_BACKEND,
                pisa_id=pisa_id,
                chains=chains,
                pisa_chain_ids=chains,
                interface_area_A2=interface_area,
                solvation_energy_kcal_mol=solvation_energy,
                stabilization_energy_kcal_mol=None,
                pvalue=pvalue,
                residue_ids=tuple(dict.fromkeys(residue_ids)),
                residue_expressions=tuple(dict.fromkeys(residue_expressions)),
                annotations=annotations,
            )
        )
    return tuple(records)


def _interface_residue_entries(element: ET.Element) -> tuple[tuple[str, str, str], ...]:
    entries_with_bsa: list[tuple[str, str, str]] = []
    fallback_entries: list[tuple[str, str, str]] = []
    for molecule in element.findall("molecule"):
        chain = _normalize_pisa_chain_id(_text(molecule, "chain_id"))
        if chain is None:
            continue
        for residue_element in molecule.findall(".//residue"):
            resid = _text(residue_element, "seq_num")
            if resid is None:
                continue
            ins_code = _text(residue_element, "ins_code")
            if ins_code:
                resid = f"{resid}{ins_code}"
            resname = (_text(residue_element, "name") or "_").upper()
            entry = (chain, resid, resname)
            fallback_entries.append(entry)
            bsa = _float_text(residue_element, "bsa")
            if bsa is not None and bsa > 0.0:
                entries_with_bsa.append(entry)
    entries = entries_with_bsa or fallback_entries
    return tuple(dict.fromkeys(entries))


def _legacy_interface_residue_entries(element: ET.Element) -> tuple[tuple[str, str, str], ...]:
    entries_with_bsa: list[tuple[str, str, str]] = []
    fallback_entries: list[tuple[str, str, str]] = []
    for residue_element in element.findall(".//RESIDUES/RESIDUE1/RESIDUE"):
        _collect_legacy_residue_entry(residue_element, entries_with_bsa, fallback_entries)
    for residue_element in element.findall(".//RESIDUES/RESIDUE2/RESIDUE"):
        _collect_legacy_residue_entry(residue_element, entries_with_bsa, fallback_entries)
    entries = entries_with_bsa or fallback_entries
    return tuple(dict.fromkeys(entries))


def _collect_legacy_residue_entry(
    residue_element: ET.Element,
    entries_with_bsa: list[tuple[str, str, str]],
    fallback_entries: list[tuple[str, str, str]],
) -> None:
    structure_text = _text(residue_element, "STRUCTURE")
    entry = _legacy_residue_entry(structure_text)
    if entry is None:
        return
    fallback_entries.append(entry)
    bsa = _float_text(residue_element, "BURIEDSURFACEAREA")
    if bsa is not None and bsa > 0.0:
        entries_with_bsa.append(entry)


def _legacy_residue_entry(text: str | None) -> tuple[str, str, str] | None:
    if text is None:
        return None
    match = re.match(
        r"\s*(?P<chain>[^:\s]+):(?P<resname>[A-Za-z0-9]{1,3})\s*(?P<resid>-?\d+[A-Za-z]?)",
        text,
    )
    if match is None:
        return None
    return (
        _normalize_pisa_chain_id(match.group("chain")) or "_",
        match.group("resid"),
        match.group("resname").upper(),
    )


def _hydrate_biological_interface_store(
    store: BiologicalInterfaceStore,
    residues_by_key: dict[tuple[str, str], Any],
) -> BiologicalInterfaceStore:
    records: list[BiologicalInterfaceRecord] = []
    for record in store.records:
        residue_ids: list[str] = []
        expressions: list[str] = []
        chains: list[str] = []
        for expression in record.residue_expressions:
            chain, resid = _chain_resid_from_expression(expression)
            residue = residues_by_key.get((chain, resid))
            if residue is None:
                continue
            residue_ids.append(object_id_for_residue(residue))
            expressions.append(_residue_expression(residue))
            chains.append(chain_for_residue(residue))
        unique_residue_ids = tuple(dict.fromkeys(residue_ids))
        unique_expressions = tuple(dict.fromkeys(expressions))
        if not unique_expressions:
            continue
        unique_chains = tuple(dict.fromkeys(chains)) or record.chains
        records.append(
            BiologicalInterfaceRecord(
                id=record.id,
                name=record.name,
                backend=record.backend,
                pisa_id=record.pisa_id,
                chains=unique_chains,
                pisa_chain_ids=record.pisa_chain_ids,
                interface_area_A2=record.interface_area_A2,
                solvation_energy_kcal_mol=record.solvation_energy_kcal_mol,
                stabilization_energy_kcal_mol=record.stabilization_energy_kcal_mol,
                pvalue=record.pvalue,
                residue_ids=unique_residue_ids,
                residue_expressions=unique_expressions,
                annotations={
                    **record.annotations,
                    "chains": list(unique_chains),
                    "participating_residues": list(unique_residue_ids),
                    "participating_residue_count": len(unique_residue_ids),
                },
            )
        )
    return BiologicalInterfaceStore(
        backend=store.backend,
        records=tuple(records),
        limitations=store.limitations,
    )


def _record_with_input_provenance(
    record: BiologicalInterfaceRecord,
    provenance: str,
) -> BiologicalInterfaceRecord:
    return BiologicalInterfaceRecord(
        id=record.id,
        name=record.name,
        backend=record.backend,
        pisa_id=record.pisa_id,
        chains=record.chains,
        pisa_chain_ids=record.pisa_chain_ids,
        interface_area_A2=record.interface_area_A2,
        solvation_energy_kcal_mol=record.solvation_energy_kcal_mol,
        stabilization_energy_kcal_mol=record.stabilization_energy_kcal_mol,
        pvalue=record.pvalue,
        residue_ids=record.residue_ids,
        residue_expressions=record.residue_expressions,
        annotations={**record.annotations, "input_provenance": provenance},
    )


def _pisa_config_args(temp_path: Path) -> tuple[list[str], tuple[str, ...]]:
    env_config = os.environ.get("MOLINSPECT_PISA_CONFIG") or os.environ.get("PISA_CONF_FILE")
    if env_config:
        config = Path(env_config).expanduser()
        if config.exists():
            return [str(config)], ()
        return [], (f"PISA config path does not exist: {config}",)

    setup_dir = os.environ.get("PISA_SETUP_DIR")
    if not setup_dir:
        return [], ()
    setup_path = Path(setup_dir).expanduser()
    template = setup_path / "pisa_cfg_tmp"
    if not template.exists():
        return [], (f"PISA_SETUP_DIR does not contain pisa_cfg_tmp: {setup_path}",)
    output = temp_path / "pisa.cfg"
    text = template.read_text()
    text = text.replace("path_dataroot", str(temp_path)).replace("path_to_setup", str(setup_path))
    output.write_text(text)
    return [str(output)], ()


def _pisa_command() -> str | None:
    env_command = os.environ.get("MOLINSPECT_PISA_COMMAND")
    if env_command:
        return env_command
    return shutil.which("pisa")


def _residues_by_chain_resid(universe: Any) -> dict[tuple[str, str], Any]:
    return {
        (chain_for_residue(residue), _residue_number(residue)): residue
        for residue in universe.residues
    }


def _residue_expression(residue: Any) -> str:
    chain = chain_for_residue(residue)
    resid = _residue_number(residue)
    resname = str(getattr(residue, "resname", "")).strip()
    return _residue_expression_from_parts(chain, resid, resname)


def _residue_expression_from_parts(chain: str, resid: str, resname: str | None) -> str:
    parts = [f"chain {chain}", f"resid {resid}"]
    if resname and resname != "_":
        parts.append(f"resname {resname}")
    return " and ".join(parts)


def _chain_resid_from_expression(expression: str) -> tuple[str, str]:
    parts = expression.split()
    try:
        return parts[1], parts[4]
    except IndexError:
        return "_", "_"


def _residue_number(residue: Any) -> str:
    resid = str(getattr(residue, "resid", "")).strip() or "_"
    icode = str(getattr(residue, "icode", "")).strip()
    return f"{resid}{icode}" if icode else resid


def _normalize_pisa_chain_id(raw_chain_id: str | None) -> str | None:
    if raw_chain_id is None:
        return None
    cleaned = raw_chain_id.strip()
    if not cleaned:
        return None
    ligand_match = re.search(r"\]([^:\-\s]+)", cleaned)
    if ligand_match is not None:
        return ligand_match.group(1) or None
    return cleaned.split(":", 1)[0].split("-", 1)[0] or None


def _pisa_session_name(source: Path) -> str:
    cleaned = re.sub(r"[^A-Za-z0-9]", "", source.stem)
    return f"molinspect{cleaned[:20] or 'target'}"


def _bond_counts(element: ET.Element) -> dict[str, int]:
    counts: dict[str, int] = {}
    for xml_tag, label in (
        ("h-bonds", "hydrogen_bonds"),
        ("salt-bridges", "salt_bridges"),
        ("ss-bonds", "disulfide_bonds"),
        ("cov-bonds", "covalent_bonds"),
    ):
        parent = element.find(xml_tag)
        if parent is None:
            continue
        count = _int_text(parent, "n_bonds")
        if count is not None:
            counts[label] = count
    return counts


def _legacy_bond_counts(element: ET.Element) -> dict[str, int]:
    counts: dict[str, int] = {}
    for xml_tag, label in (
        ("HYDROGENBONDS", "hydrogen_bonds"),
        ("SALTBRIDGES", "salt_bridges"),
        ("DISULFIDEBONDS", "disulfide_bonds"),
        ("COVALENTBONDS", "covalent_bonds"),
    ):
        parent = element.find(xml_tag)
        if parent is None:
            continue
        counts[label] = len(parent.findall("STRUCTURE"))
    return {key: value for key, value in counts.items() if value}


def _legacy_interface_area(element: ET.Element) -> float | None:
    areas = _legacy_partner_areas(element)
    if not areas:
        return None
    return round(sum(areas) / len(areas), 3)


def _legacy_partner_areas(element: ET.Element) -> list[float]:
    areas: list[float] = []
    for tag in ("INTERFACESUMMARY/STRUCTURE1/SOLVENTAREA1", "INTERFACESUMMARY/STRUCTURE2/SOLVENTAREA2"):
        area = _float_text(element, f"{tag}/INTERFACEAREA")
        if area is not None:
            areas.append(area)
    return areas


def _legacy_solvation_energy(element: ET.Element) -> float | None:
    values: list[float] = []
    for tag in (
        "INTERFACESUMMARY/STRUCTURE1/SOLVATIONENERGY1",
        "INTERFACESUMMARY/STRUCTURE2/SOLVATIONENERGY2",
    ):
        value = _float_text(element, f"{tag}/GAINCOMPLEXFORMATION")
        if value is not None:
            values.append(value)
    if not values:
        return None
    return round(sum(values), 3)


def _legacy_pvalue(element: ET.Element) -> float | None:
    values: list[float] = []
    for tag in (
        "INTERFACESUMMARY/STRUCTURE1/SOLVATIONENERGY1",
        "INTERFACESUMMARY/STRUCTURE2/SOLVATIONENERGY2",
    ):
        value = _float_text(element, f"{tag}/PVALUE")
        if value is not None:
            values.append(value)
    if not values:
        return None
    return round(max(values), 3)


def _text(element: ET.Element, tag: str) -> str | None:
    child = element.find(tag)
    if child is None or child.text is None:
        return None
    text = child.text.strip()
    return text or None


def _float_text(element: ET.Element, tag: str) -> float | None:
    text = _text(element, tag)
    if text is None:
        return None
    try:
        return float(text)
    except ValueError:
        return None


def _int_text(element: ET.Element, tag: str) -> int | None:
    text = _text(element, tag)
    if text is None:
        return None
    try:
        return int(text)
    except ValueError:
        return None


def _id_part(value: Any) -> str:
    cleaned = str(value).strip() or "_"
    return cleaned.replace(":", "_").replace(" ", "_")
