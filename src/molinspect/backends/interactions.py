"""Optional interaction backends with auditable fallback behavior."""

from __future__ import annotations

import json
import shutil
import subprocess
import sys
import tempfile
import xml.etree.ElementTree as ET
from dataclasses import dataclass
from functools import lru_cache
from pathlib import Path
from typing import Any

from .common import (
    module_can_import,
    short_process_error,
    source_cache_key,
    static_structure_source,
)
from ..metrics import closest_heavy_atom_pair
from ..relations import chain_for_residue, relation_priority_for_type
from ..schemas import EvidenceItem, Relation
from ..definitions import (
    CONTACT_CUTOFF_A,
    HBOND_DISTANCE_CUTOFF_A,
    HYDROPHOBIC_CUTOFF_A,
    METAL_COORDINATION_CUTOFF_A,
    MOLINSPECT_HEURISTIC_BACKEND,
    PI_STACKING_CENTER_CUTOFF_A,
    SALT_BRIDGE_CUTOFF_A,
    STERIC_CLASH_MIN_VDW_OVERLAP_A,
    definition,
)

ARPEGGIO_BACKEND = "PDBe Arpeggio"
ARPEGGIO_TIMEOUT_S = 8
PLIP_BACKEND = "PLIP"
PLIP_TIMEOUT_S = 60

INTERACTION_BACKEND_FOCUS = frozenset(
    {
        "contacts",
        "metal_coordination",
        "hydrogen_bonds",
        "salt_bridges",
        "hydrophobic_contacts",
        "pi_stacking",
        "steric_clashes",
    }
)

_ARPEGGIO_CONTACT_TO_DEFINITION = {
    "metal_complex": "metal_coordination",
    "ionic": "salt_bridge",
    "hbond": "hydrogen_bond",
    "weak_hbond": "polar_contact_candidate",
    "polar": "polar_contact_candidate",
    "aromatic": "pi_stacking",
    "hydrophobic": "hydrophobic_contact",
    "clash": "steric_clash",
    "vdw_clash": "steric_clash",
    "vdw": "nonbonded_contact",
    "proximal": "near",
}

_ARPEGGIO_CONTACT_PRIORITY = (
    "metal_complex",
    "ionic",
    "hbond",
    "aromatic",
    "hydrophobic",
    "weak_hbond",
    "polar",
    "clash",
    "vdw_clash",
    "vdw",
    "proximal",
)


@dataclass(frozen=True, slots=True)
class AtomKey:
    """Stable atom identity shared by MDAnalysis and backend JSON outputs."""

    chain: str
    resid: str
    icode: str
    resname: str
    atom_name: str

    @property
    def label(self) -> str:
        resid = f"{self.resid}{self.icode}" if self.icode else self.resid
        return f"{self.chain}:{resid}:{self.resname}:{self.atom_name}"


@dataclass(frozen=True, slots=True)
class BackendContact:
    """One backend contact mapped onto MolInspect's relation vocabulary."""

    bgn_keys: tuple[AtomKey, ...]
    end_keys: tuple[AtomKey, ...]
    definition_id: str
    backend_contacts: tuple[str, ...]
    distance_A: float
    contact_type: str
    interacting_entities: str | None
    backend: str
    input_provenance: str | None = None
    backend_priority: int = 1


@dataclass(frozen=True, slots=True)
class InteractionStore:
    """Interaction contacts indexed for selected-object lookups."""

    backend: str
    contacts: tuple[BackendContact, ...] = ()
    limitations: tuple[str, ...] = ()
    input_provenance: str | None = None

    def relation_between(
        self,
        source_atoms: Any,
        target_atoms: Any,
        source_id: str,
        target_id: str,
    ) -> Relation | None:
        """Return the strongest backend contact between two atom groups."""

        if not self.contacts:
            return None

        source_keys = {_atom_key_for_mda_atom(atom) for atom in source_atoms}
        target_keys = {_atom_key_for_mda_atom(atom) for atom in target_atoms}
        candidates: list[tuple[int, int, float, BackendContact, bool]] = []
        for contact in self.contacts:
            bgn_in_source = bool(source_keys.intersection(contact.bgn_keys))
            end_in_source = bool(source_keys.intersection(contact.end_keys))
            bgn_in_target = bool(target_keys.intersection(contact.bgn_keys))
            end_in_target = bool(target_keys.intersection(contact.end_keys))
            if bgn_in_source and end_in_target:
                candidates.append(
                    (
                        relation_priority_for_type(contact.definition_id),
                        contact.backend_priority,
                        contact.distance_A,
                        contact,
                        True,
                    )
                )
            elif end_in_source and bgn_in_target:
                candidates.append(
                    (
                        relation_priority_for_type(contact.definition_id),
                        contact.backend_priority,
                        contact.distance_A,
                        contact,
                        False,
                    )
                )

        if not candidates:
            return None

        _, _, _, best_contact, source_is_bgn = sorted(candidates, key=lambda item: item[:3])[0]
        return _relation_from_backend_contact(best_contact, source_id, target_id, source_is_bgn, self)


def build_interaction_store(source_files: tuple[Path, ...]) -> InteractionStore:
    """Build the best available interaction store for a static source."""

    source = static_structure_source(source_files)
    if source is None:
        return InteractionStore(
            backend=ARPEGGIO_BACKEND,
            limitations=("PDBe Arpeggio interactions require a static PDB/mmCIF source.",),
        )
    source_key = source_cache_key(source)
    stores = [_cached_plip_store(source_key), _cached_arpeggio_store(source_key)]
    contacts = tuple(contact for store in stores for contact in store.contacts)
    limitations = tuple(
        limitation for store in stores for limitation in store.limitations if limitation
    )
    backends = tuple(dict.fromkeys(contact.backend for contact in contacts))
    if contacts:
        return InteractionStore(
            backend="+".join(backends),
            contacts=contacts,
            limitations=limitations,
        )
    return InteractionStore(
        backend=ARPEGGIO_BACKEND,
        limitations=limitations or ("No interaction backend produced usable contacts.",),
    )


def interaction_backend_is_available() -> bool:
    """Return whether PDBe Arpeggio can run in this environment."""

    return (
        shutil.which("pdbe-arpeggio") is not None
        and module_can_import("arpeggio")
        and module_can_import("openbabel")
        and module_can_import("gemmi")
    )


def plip_is_available() -> bool:
    """Return whether PLIP is importable or exposed as a command."""

    return module_can_import("plip") or shutil.which("plip") is not None


@lru_cache(maxsize=16)
def _cached_arpeggio_store(source_key: tuple[str, int, int]) -> InteractionStore:
    source = Path(source_key[0])
    if not interaction_backend_is_available():
        return InteractionStore(
            backend=ARPEGGIO_BACKEND,
            limitations=(
                "PDBe Arpeggio is unavailable; install pdbe-arpeggio plus the openbabel Python "
                "module to enable backend-typed interactions.",
            ),
        )

    with tempfile.TemporaryDirectory(prefix="molinspect-arpeggio-") as temp_dir:
        temp_path = Path(temp_dir)
        try:
            arpeggio_input, input_provenance = _arpeggio_input(source, temp_path)
        except Exception as exc:
            return InteractionStore(
                backend=ARPEGGIO_BACKEND,
                limitations=(f"PDBe Arpeggio input preparation failed: {exc}",),
            )

        output_dir = temp_path / "out"
        output_dir.mkdir()
        command = [
            sys.executable,
            "-m",
            "arpeggio.scripts.process_protein_cli",
            str(arpeggio_input),
            "-o",
            str(output_dir),
            "-m",
        ]
        try:
            completed = subprocess.run(
                command,
                check=False,
                capture_output=True,
                text=True,
                timeout=ARPEGGIO_TIMEOUT_S,
            )
        except subprocess.TimeoutExpired:
            return InteractionStore(
                backend=ARPEGGIO_BACKEND,
                limitations=(f"PDBe Arpeggio timed out after {ARPEGGIO_TIMEOUT_S} seconds.",),
            )

        if completed.returncode != 0:
            detail = short_process_error(completed.stderr or completed.stdout)
            return InteractionStore(
                backend=ARPEGGIO_BACKEND,
                limitations=(f"PDBe Arpeggio failed: {detail}",),
            )

        json_paths = sorted(output_dir.glob("*.json"))
        if not json_paths:
            return InteractionStore(
                backend=ARPEGGIO_BACKEND,
                limitations=("PDBe Arpeggio completed but produced no contact JSON.",),
            )
        try:
            raw_contacts = json.loads(json_paths[0].read_text())
        except Exception as exc:
            return InteractionStore(
                backend=ARPEGGIO_BACKEND,
                limitations=(f"PDBe Arpeggio contact JSON could not be parsed: {exc}",),
            )

    contacts = tuple(_iter_arpeggio_contacts(raw_contacts))
    if not contacts:
        return InteractionStore(
            backend=ARPEGGIO_BACKEND,
            limitations=("PDBe Arpeggio produced no MolInspect-mapped interaction contacts.",),
            input_provenance=input_provenance,
        )
    return InteractionStore(
        backend=ARPEGGIO_BACKEND,
        contacts=contacts,
        input_provenance=input_provenance,
    )


def _iter_arpeggio_contacts(raw_contacts: Any) -> list[BackendContact]:
    contacts: list[BackendContact] = []
    if not isinstance(raw_contacts, list):
        return contacts

    for entry in raw_contacts:
        if not isinstance(entry, dict):
            continue
        contact_names = _contact_names(entry.get("contact"))
        definition_id = _definition_id_for_arpeggio_contacts(contact_names)
        if definition_id is None:
            continue
        bgn_keys = _endpoint_keys(entry.get("bgn"))
        end_keys = _endpoint_keys(entry.get("end"))
        if not bgn_keys or not end_keys:
            continue
        try:
            distance_A = round(float(entry["distance"]), 3)
        except (KeyError, TypeError, ValueError):
            continue
        contacts.append(
            BackendContact(
                bgn_keys=bgn_keys,
                end_keys=end_keys,
                definition_id=definition_id,
                backend_contacts=contact_names,
                distance_A=distance_A,
                contact_type=str(entry.get("type") or ""),
                interacting_entities=str(entry.get("interacting_entities") or "")
                or None,
                backend=ARPEGGIO_BACKEND,
                input_provenance="PDBe Arpeggio; source mmCIF or Gemmi-converted PDB",
                backend_priority=1,
            )
        )
    return contacts


@lru_cache(maxsize=16)
def _cached_plip_store(source_key: tuple[str, int, int]) -> InteractionStore:
    source = Path(source_key[0])
    if not plip_is_available() or not _source_has_nonwater_hetatm(source):
        return InteractionStore(backend=PLIP_BACKEND)

    with tempfile.TemporaryDirectory(prefix="molinspect-plip-") as temp_dir:
        output_dir = Path(temp_dir) / "out"
        output_dir.mkdir()
        command = [
            shutil.which("plip") or "plip",
            "-f",
            str(source),
            "-o",
            str(output_dir),
            "-x",
            "-q",
            "--maxthreads",
            "1",
            "--name",
            "plip",
        ]
        try:
            completed = subprocess.run(
                command,
                check=False,
                capture_output=True,
                text=True,
                timeout=PLIP_TIMEOUT_S,
            )
        except subprocess.TimeoutExpired:
            return InteractionStore(
                backend=PLIP_BACKEND,
                limitations=(f"PLIP timed out after {PLIP_TIMEOUT_S} seconds.",),
            )
        if completed.returncode != 0:
            return InteractionStore(
                backend=PLIP_BACKEND,
                limitations=(f"PLIP failed: {short_process_error(completed.stderr)}",),
            )
        xml_paths = sorted(output_dir.glob("*.xml"))
        if not xml_paths:
            return InteractionStore(
                backend=PLIP_BACKEND,
                limitations=("PLIP completed but produced no XML report.",),
            )
        atom_keys = _pdb_atom_keys_by_serial(source)
        contacts = tuple(_iter_plip_contacts(xml_paths[0], atom_keys))
        return InteractionStore(backend=PLIP_BACKEND, contacts=contacts)


def _iter_plip_contacts(
    xml_path: Path,
    atom_keys_by_serial: dict[str, AtomKey],
) -> list[BackendContact]:
    root = ET.parse(xml_path).getroot()
    contacts: list[BackendContact] = []
    for tag, definition_id, source_fields, target_fields, distance_field in (
        (
            "hydrophobic_interaction",
            "hydrophobic_contact",
            ("ligcarbonidx",),
            ("protcarbonidx",),
            "dist",
        ),
        ("hydrogen_bond", "hydrogen_bond", ("donor_idx",), ("acceptor_idx",), "dist_da"),
        ("salt_bridge", "salt_bridge", ("lig_idx_list",), ("prot_idx_list",), "dist"),
        ("pi_stack", "pi_stacking", ("lig_idx_list",), ("prot_idx_list",), "centdist"),
        ("metal_complex", "metal_coordination", ("metal_idx",), ("target_idx",), "dist"),
    ):
        for element in root.findall(f".//{tag}"):
            bgn_keys = _plip_atom_keys(element, source_fields, atom_keys_by_serial)
            end_keys = _plip_atom_keys(element, target_fields, atom_keys_by_serial)
            distance_A = _float_text(element, distance_field)
            if not bgn_keys or not end_keys or distance_A is None:
                continue
            contacts.append(
                BackendContact(
                    bgn_keys=bgn_keys,
                    end_keys=end_keys,
                    definition_id=definition_id,
                    backend_contacts=(tag,),
                    distance_A=round(distance_A, 3),
                    contact_type=tag,
                    interacting_entities="protein_ligand",
                    backend=PLIP_BACKEND,
                    input_provenance="PLIP XML report",
                    backend_priority=0,
                )
            )
    return contacts


def _plip_atom_keys(
    element: ET.Element,
    fields: tuple[str, ...],
    atom_keys_by_serial: dict[str, AtomKey],
) -> tuple[AtomKey, ...]:
    keys: list[AtomKey] = []
    for field in fields:
        child = element.find(field)
        if child is None:
            continue
        if field.endswith("_list"):
            for index_element in child.findall("idx"):
                key = atom_keys_by_serial.get((index_element.text or "").strip())
                if key is not None:
                    keys.append(key)
        else:
            key = atom_keys_by_serial.get((child.text or "").strip())
            if key is not None:
                keys.append(key)
    return tuple(dict.fromkeys(keys))


def _pdb_atom_keys_by_serial(source: Path) -> dict[str, AtomKey]:
    atom_keys: dict[str, AtomKey] = {}
    if source.suffix.lower() not in {".pdb", ".ent"}:
        return atom_keys
    with source.open() as handle:
        for line in handle:
            if not line.startswith(("ATOM", "HETATM")):
                continue
            serial = line[6:11].strip()
            if not serial:
                continue
            chain = line[21].strip() or "_"
            resid = line[22:26].strip()
            icode = line[26].strip()
            resname = line[17:20].strip().upper()
            atom_name = line[12:16].strip()
            atom_keys[serial] = AtomKey(
                chain=chain,
                resid=resid,
                icode=icode,
                resname=resname,
                atom_name=atom_name,
            )
    return atom_keys


def _source_has_nonwater_hetatm(source: Path) -> bool:
    if source.suffix.lower() not in {".pdb", ".ent"}:
        return True
    waters = {"HOH", "WAT", "H2O", "TIP3", "TIP3P", "SOL", "SOLV", "SPC"}
    with source.open() as handle:
        for line in handle:
            if line.startswith("HETATM") and line[17:20].strip().upper() not in waters:
                return True
    return False


def _float_text(element: ET.Element, tag: str) -> float | None:
    child = element.find(tag)
    if child is None or child.text is None:
        return None
    try:
        return float(child.text)
    except ValueError:
        return None


def _relation_from_backend_contact(
    contact: BackendContact,
    source_id: str,
    target_id: str,
    source_is_bgn: bool,
    store: InteractionStore,
) -> Relation:
    definition_record = definition(contact.definition_id)
    backend_reference = "plip_docs" if contact.backend == PLIP_BACKEND else "arpeggio_paper"
    source_keys = contact.bgn_keys if source_is_bgn else contact.end_keys
    target_keys = contact.end_keys if source_is_bgn else contact.bgn_keys
    evidence = [
        EvidenceItem(
            type="metric",
            metric="backend_contact_distance",
            value=contact.distance_A,
            unit="angstrom",
            source=f"{source_keys[0].label}->{target_keys[0].label}",
        ),
        EvidenceItem(
            type="method",
            metric="interaction_backend",
            value=contact.backend,
            source=contact.input_provenance or store.input_provenance,
        ),
        EvidenceItem(
            type="method",
            metric="backend_contact_types",
            value=list(contact.backend_contacts),
            source=contact.interacting_entities,
        ),
    ]
    return Relation(
        source=source_id,
        target=target_id,
        type=contact.definition_id,
        category=definition_record.category,
        confidence="backend",
        backend=contact.backend,
        definition_id=contact.definition_id,
        definition_source="backend",
        reference_keys=list(dict.fromkeys((*definition_record.reference_keys, backend_reference))),
        min_distance_A=contact.distance_A,
        cutoff_A=_cutoff_for_relation(contact.definition_id),
        source_atom=source_keys[0].label,
        target_atom=target_keys[0].label,
        evidence=evidence,
        limitations=list(definition_record.limitations),
    )


def _arpeggio_input(source: Path, temp_dir: Path) -> tuple[Path, str]:
    suffix = source.suffix.lower()
    if suffix in {".cif", ".mmcif"}:
        return source, "source mmCIF"
    if suffix != ".pdb":
        raise ValueError(f"unsupported static structure format {suffix!r}")

    import gemmi  # type: ignore[import-not-found]

    converted = temp_dir / f"{source.stem}.cif"
    structure = gemmi.read_structure(str(source))
    structure.make_mmcif_document().write_file(str(converted))
    return converted, "Gemmi-converted mmCIF from PDB source"


def _endpoint_keys(endpoint: Any) -> tuple[AtomKey, ...]:
    if not isinstance(endpoint, dict):
        return ()
    atom_field = str(endpoint.get("auth_atom_id") or "").strip()
    atom_names = [part.strip() for part in atom_field.split(",") if part.strip()]
    if not atom_names:
        return ()
    chain = str(endpoint.get("auth_asym_id") or "_").strip() or "_"
    resid = str(endpoint.get("auth_seq_id") or "").strip()
    resname = str(endpoint.get("label_comp_id") or "").strip().upper()
    icode = str(endpoint.get("pdbx_PDB_ins_code") or "").strip()
    if icode in {".", "?"}:
        icode = ""
    return tuple(
        AtomKey(
            chain=chain,
            resid=resid,
            icode=icode,
            resname=resname,
            atom_name=atom_name,
        )
        for atom_name in atom_names
    )


def _atom_key_for_mda_atom(atom: Any) -> AtomKey:
    residue = atom.residue
    icode = str(getattr(residue, "icode", "") or "").strip()
    if icode in {".", "?"}:
        icode = ""
    return AtomKey(
        chain=chain_for_residue(residue),
        resid=str(getattr(residue, "resid", "")).strip(),
        icode=icode,
        resname=str(getattr(residue, "resname", "")).strip().upper(),
        atom_name=str(getattr(atom, "name", "")).strip(),
    )


def _contact_names(raw_contact: Any) -> tuple[str, ...]:
    if isinstance(raw_contact, str):
        values = [raw_contact]
    elif isinstance(raw_contact, list | tuple | set):
        values = [str(value) for value in raw_contact]
    else:
        values = []
    return tuple(value.strip().lower() for value in values if value and str(value).strip())


def _definition_id_for_arpeggio_contacts(contact_names: tuple[str, ...]) -> str | None:
    contact_set = set(contact_names)
    for contact_name in _ARPEGGIO_CONTACT_PRIORITY:
        if contact_name in contact_set:
            return _ARPEGGIO_CONTACT_TO_DEFINITION[contact_name]
    return None


def _cutoff_for_relation(definition_id: str) -> float:
    if definition_id == "metal_coordination":
        return METAL_COORDINATION_CUTOFF_A
    if definition_id == "salt_bridge":
        return SALT_BRIDGE_CUTOFF_A
    if definition_id in {"hydrogen_bond", "polar_contact_candidate"}:
        return HBOND_DISTANCE_CUTOFF_A
    if definition_id == "pi_stacking":
        return PI_STACKING_CENTER_CUTOFF_A
    if definition_id == "hydrophobic_contact":
        return HYDROPHOBIC_CUTOFF_A
    if definition_id == "steric_clash":
        return STERIC_CLASH_MIN_VDW_OVERLAP_A
    return CONTACT_CUTOFF_A


def backend_relation_or_fallback(
    store: InteractionStore,
    source_atoms: Any,
    target_atoms: Any,
    source_id: str,
    target_id: str,
    fallback: Relation,
) -> Relation:
    """Return backend relation when available, otherwise the supplied fallback."""

    if store.backend == MOLINSPECT_HEURISTIC_BACKEND:
        return fallback
    relation = store.relation_between(source_atoms, target_atoms, source_id, target_id)
    if relation is not None:
        return relation
    if fallback.min_distance_A is None:
        pair = closest_heavy_atom_pair(source_atoms, target_atoms)
        if pair is not None:
            return fallback.model_copy(update={"min_distance_A": pair.distance_A})
    return fallback
