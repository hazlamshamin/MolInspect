"""Backend adapters for residue-level annotations."""

from __future__ import annotations

import shutil
import subprocess
import tempfile
from dataclasses import dataclass
from functools import lru_cache
from pathlib import Path
from typing import Any

from .common import static_structure_source
from ..definitions import (
    FREESASA_BURIED_RELATIVE_SASA_MAX,
    FREESASA_EXPOSED_RELATIVE_SASA_MIN,
    SECONDARY_STRUCTURE_NAMES,
)


@dataclass(frozen=True, slots=True)
class SecondaryStructureAnnotation:
    """One DSSP secondary-structure assignment keyed by chain/residue."""

    code: str
    name: str


@dataclass(frozen=True, slots=True)
class SasaProfile:
    """One FreeSASA residue exposure profile."""

    sasa_A2: float
    relative_sasa: float | None
    exposure: str | None


def freesasa_profiles(source_files: tuple[Path, ...]) -> tuple[dict[tuple[str, str], SasaProfile], list[str]]:
    """Return FreeSASA residue profiles plus limitations."""

    source = static_structure_source(source_files)
    if source is None:
        return {}, ["FreeSASA exact SASA was not run because no PDB/mmCIF source was available."]

    try:
        import freesasa  # type: ignore[import-not-found]
    except ImportError:
        return {}, ["FreeSASA exact SASA is unavailable because freesasa is not installed."]

    try:
        structure = freesasa.Structure(str(source))
        result = freesasa.calc(structure)
        residue_areas = result.residueAreas()
    except Exception as exc:
        return {}, [f"FreeSASA exact SASA could not be computed: {exc}"]

    profiles = parse_freesasa_residue_areas(residue_areas)
    if not profiles:
        return {}, ["FreeSASA exact SASA produced no residue-level areas."]
    return profiles, []


def parse_freesasa_residue_areas(residue_areas: Any) -> dict[tuple[str, str], SasaProfile]:
    """Parse FreeSASA residueAreas() output into normalized profiles."""

    profiles: dict[tuple[str, str], SasaProfile] = {}
    for raw_chain, residues in _iter_mapping_items(residue_areas):
        chain = str(raw_chain).strip() or "_"
        for raw_resid, area in _iter_mapping_items(residues):
            resid = str(raw_resid).strip()
            if not resid:
                continue
            total = _area_value(area, ("total", "totalArea", "sasa_A2"))
            if total is None:
                continue
            relative = _normalize_relative_sasa(
                _area_value(area, ("relativeTotal", "relative_total", "relative_sasa"))
            )
            profiles[(chain, resid)] = SasaProfile(
                sasa_A2=round(total, 3),
                relative_sasa=round(relative, 3) if relative is not None else None,
                exposure=_exposure_from_relative_sasa(relative),
            )
    return profiles


def dssp_annotations(
    source_files: tuple[Path, ...],
) -> tuple[dict[tuple[str, str], SecondaryStructureAnnotation], list[str]]:
    """Return DSSP residue assignments plus limitations."""

    source = static_structure_source(source_files)
    if source is None:
        return {}, ["DSSP secondary structure was not run because no PDB/mmCIF source was available."]
    annotations, limitations = _cached_dssp_annotations(str(source.resolve()))
    return annotations, list(limitations)


@lru_cache(maxsize=16)
def _cached_dssp_annotations(
    source_path: str,
) -> tuple[dict[tuple[str, str], SecondaryStructureAnnotation], tuple[str, ...]]:
    executable = shutil.which("mkdssp")
    if executable is None:
        return {}, ("DSSP secondary structure is unavailable because mkdssp is not installed.",)

    try:
        with tempfile.NamedTemporaryFile(suffix=".dssp") as output:
            subprocess.run(
                [executable, "--output-format", "dssp", source_path, output.name],
                check=True,
                capture_output=True,
                text=True,
            )
            text = Path(output.name).read_text()
    except Exception as exc:
        return {}, (f"DSSP secondary structure could not be computed: {exc}",)

    return parse_dssp(text), ()


def parse_dssp(text: str) -> dict[tuple[str, str], SecondaryStructureAnnotation]:
    """Parse DSSP fixed-width output into residue assignments."""

    annotations: dict[tuple[str, str], SecondaryStructureAnnotation] = {}
    in_table = False
    for line in text.splitlines():
        if line.startswith("  #  RESIDUE"):
            in_table = True
            continue
        if not in_table or len(line) < 17:
            continue
        resid = line[5:10].strip()
        insertion_code = line[10].strip()
        chain = line[11].strip() or "_"
        code = line[16] if line[16].strip() else " "
        if not resid or line[13] == "!":
            continue
        residue_id = f"{resid}{insertion_code}" if insertion_code else resid
        annotations[(chain, residue_id)] = SecondaryStructureAnnotation(
            code=code,
            name=SECONDARY_STRUCTURE_NAMES.get(code, "other"),
        )
    return annotations


def dssp_key_for_residue(residue: Any) -> tuple[str, str]:
    """Return the chain/residue key used by DSSP and FreeSASA maps."""

    chain = _chain_for_residue(residue)
    resid = str(getattr(residue, "resid", "")).strip()
    insertion_code = str(getattr(residue, "icode", "")).strip()
    return chain, f"{resid}{insertion_code}" if insertion_code else resid


def _iter_mapping_items(value: Any) -> Any:
    if hasattr(value, "items"):
        return value.items()
    return ()


def _area_value(area: Any, names: tuple[str, ...]) -> float | None:
    for name in names:
        raw_value = getattr(area, name, None)
        if raw_value is None and isinstance(area, dict):
            raw_value = area.get(name)
        if callable(raw_value):
            raw_value = raw_value()
        if raw_value is None:
            continue
        try:
            return float(raw_value)
        except (TypeError, ValueError):
            continue
    return None


def _normalize_relative_sasa(value: float | None) -> float | None:
    if value is None:
        return None
    if value > 1.0:
        return value / 100.0
    return value


def _exposure_from_relative_sasa(relative_sasa: float | None) -> str | None:
    if relative_sasa is None:
        return None
    if relative_sasa <= FREESASA_BURIED_RELATIVE_SASA_MAX:
        return "buried"
    if relative_sasa >= FREESASA_EXPOSED_RELATIVE_SASA_MIN:
        return "exposed"
    return "partially_buried"


def _chain_for_residue(residue: Any) -> str:
    for atom in residue.atoms:
        chain = str(getattr(atom, "chainID", "")).strip()
        if chain:
            return chain
    segment = str(getattr(residue, "segid", "")).strip()
    return segment or "_"
