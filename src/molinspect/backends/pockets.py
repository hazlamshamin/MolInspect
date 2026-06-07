"""Optional pocket detection backends.

P2Rank is the preferred backend because it predicts ligand-binding sites and
reports calibrated probabilities. fpocket remains a local geometric fallback;
its failures are surfaced as limitations rather than silently rebranded as
true pocket evidence.
"""

from __future__ import annotations

import csv
import os
import re
import shutil
import subprocess
import tempfile
from dataclasses import dataclass
from functools import lru_cache
from pathlib import Path
from typing import Any

from .common import short_process_error, source_cache_key, static_structure_source
from ..objects import object_id_for_residue
from ..relations import chain_for_residue
from ..schemas import ObjectRef
from ..definitions import definition

P2RANK_BACKEND = "P2Rank"
FPOCKET_BACKEND = "fpocket"
POCKET_TIMEOUT_S = 30

_P2RANK_RESIDUE_PATTERN = re.compile(r"(?P<chain>[^_\s]+)_(?P<resid>[+-]?\d+[A-Za-z]?)")


@dataclass(frozen=True, slots=True)
class PocketRecord:
    """Backend pocket result mapped to residue selections."""

    id: str
    name: str
    backend: str
    rank: int
    center_A: tuple[float, float, float] | None
    score: float | None
    probability: float | None
    residue_ids: tuple[str, ...]
    residue_expressions: tuple[str, ...]
    annotations: dict[str, Any]


@dataclass(frozen=True, slots=True)
class PocketStore:
    """Detected pockets and their residue memberships."""

    backend: str | None = None
    records: tuple[PocketRecord, ...] = ()
    limitations: tuple[str, ...] = ()


def detect_pockets(universe: Any, source_files: tuple[Path, ...]) -> PocketStore:
    """Run the preferred available pocket backend and return mapped pockets."""

    source = static_structure_source(source_files)
    if source is None:
        return PocketStore(limitations=("Pocket detection requires a static PDB/mmCIF source.",))

    residues_by_key = _residues_by_chain_resid(universe)
    store = _cached_pocket_detection(source_cache_key(source), tuple(sorted(residues_by_key)))
    return _hydrate_pocket_store(store, residues_by_key)


def p2rank_is_available() -> bool:
    """Return whether a P2Rank command is discoverable."""

    return _p2rank_command() is not None


def fpocket_is_available() -> bool:
    """Return whether an fpocket command is discoverable."""

    return shutil.which("fpocket") is not None


@lru_cache(maxsize=16)
def _cached_pocket_detection(
    source_key: tuple[str, int, int],
    residue_keys: tuple[tuple[str, str], ...],
) -> PocketStore:
    source = Path(source_key[0])
    p2rank_command = _p2rank_command()
    if p2rank_command is not None:
        return _run_p2rank(source, p2rank_command)

    fpocket_command = shutil.which("fpocket")
    if fpocket_command is not None:
        fpocket_store = _run_fpocket(source, fpocket_command, residue_keys)
        if fpocket_store.records:
            return fpocket_store
        return PocketStore(
            backend=FPOCKET_BACKEND,
            limitations=(
                "P2Rank is not installed; fpocket fallback did not produce usable pockets. "
                f"{' '.join(fpocket_store.limitations)}",
            ),
        )

    return PocketStore(
        limitations=(
            "No pocket backend is available. Install P2Rank (`prank`) for ligand-binding-site "
            "predictions or fpocket for geometric pocket fallback.",
        )
    )


def _run_p2rank(source: Path, command: str) -> PocketStore:
    with tempfile.TemporaryDirectory(prefix="molinspect-p2rank-") as temp_dir:
        output_dir = Path(temp_dir) / "out"
        cmd = [
            command,
            "predict",
            "-f",
            str(source),
            "-o",
            str(output_dir),
            "-visualizations",
            "0",
            "-threads",
            "1",
        ]
        try:
            completed = subprocess.run(
                cmd,
                check=False,
                capture_output=True,
                text=True,
                timeout=POCKET_TIMEOUT_S,
            )
        except subprocess.TimeoutExpired:
            return PocketStore(
                backend=P2RANK_BACKEND,
                limitations=(f"P2Rank timed out after {POCKET_TIMEOUT_S} seconds.",),
            )
        if completed.returncode != 0:
            return PocketStore(
                backend=P2RANK_BACKEND,
                limitations=(f"P2Rank failed: {short_process_error(completed.stderr)}",),
            )
        predictions = sorted(output_dir.rglob("*_predictions.csv"))
        if not predictions:
            return PocketStore(
                backend=P2RANK_BACKEND,
                limitations=("P2Rank completed but produced no *_predictions.csv file.",),
            )
        return PocketStore(
            backend=P2RANK_BACKEND,
            records=tuple(_parse_p2rank_predictions(predictions[0])),
        )


def _parse_p2rank_predictions(path: Path) -> list[PocketRecord]:
    records: list[PocketRecord] = []
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle, skipinitialspace=True)
        for index, row in enumerate(reader, start=1):
            rank = _as_int(_row_value(row, "rank")) or index
            name = _row_value(row, "name") or f"pocket{rank}"
            residue_tokens = (_row_value(row, "residue_ids") or "").split()
            residue_expressions = tuple(
                expression
                for token in residue_tokens
                if (expression := _p2rank_residue_expression(token)) is not None
            )
            residue_ids = tuple(
                f"residue:{token.replace('_', ':')}" for token in residue_tokens if token
            )
            center = _center_from_row(row)
            pocket_id = f"pocket:p2rank:{rank}"
            annotations: dict[str, Any] = {
                "definition_id": "pocket",
                "definition_source": definition("pocket").source,
                "reference_keys": list(definition("pocket").reference_keys),
                "backend": P2RANK_BACKEND,
                "method": "p2rank_ligand_binding_site_prediction",
                "rank": rank,
                "score": _as_float(_row_value(row, "score")),
                "probability": _as_float(_row_value(row, "probability")),
                "center_A": list(center) if center is not None else None,
                "lining_residue_count": len(residue_expressions),
                "limitation": definition("pocket").limitations[0],
            }
            records.append(
                PocketRecord(
                    id=pocket_id,
                    name=f"{name}_p2rank",
                    backend=P2RANK_BACKEND,
                    rank=rank,
                    center_A=center,
                    score=annotations["score"],
                    probability=annotations["probability"],
                    residue_ids=residue_ids,
                    residue_expressions=residue_expressions,
                    annotations=annotations,
                )
            )
    return records


def _run_fpocket(
    source: Path,
    command: str,
    residue_keys: tuple[tuple[str, str], ...],
) -> PocketStore:
    with tempfile.TemporaryDirectory(prefix="molinspect-fpocket-") as temp_dir:
        temp_path = Path(temp_dir)
        copied = temp_path / source.name
        shutil.copy2(source, copied)
        try:
            completed = subprocess.run(
                [command, "-f", copied.name],
                cwd=temp_path,
                check=False,
                capture_output=True,
                text=True,
                timeout=POCKET_TIMEOUT_S,
            )
        except subprocess.TimeoutExpired:
            return PocketStore(
                backend=FPOCKET_BACKEND,
                limitations=(f"fpocket timed out after {POCKET_TIMEOUT_S} seconds.",),
            )
        if completed.returncode != 0:
            return PocketStore(
                backend=FPOCKET_BACKEND,
                limitations=(f"fpocket failed: {short_process_error(completed.stderr)}",),
            )
        output_root = temp_path / f"{source.stem}_out" / "pockets"
        pocket_files = sorted(output_root.glob("pocket*_atm.pdb"))
        if not pocket_files:
            return PocketStore(
                backend=FPOCKET_BACKEND,
                limitations=("fpocket completed but produced no pocket*_atm.pdb files.",),
            )
        known_residue_keys = set(residue_keys)
        records = [
            record
            for record in (
                _parse_fpocket_file(path, rank=index, known_residue_keys=known_residue_keys)
                for index, path in enumerate(pocket_files, start=1)
            )
            if record is not None
        ]
        return PocketStore(backend=FPOCKET_BACKEND, records=tuple(records))


def _parse_fpocket_file(
    path: Path,
    rank: int,
    known_residue_keys: set[tuple[str, str]],
) -> PocketRecord | None:
    residue_tokens: set[tuple[str, str, str]] = set()
    center_values: list[tuple[float, float, float]] = []
    with path.open() as handle:
        for line in handle:
            if not line.startswith(("ATOM", "HETATM")):
                continue
            chain = line[21].strip() or "_"
            resid = line[22:27].strip()
            resname = line[17:20].strip().upper()
            if not resid or (chain, resid) not in known_residue_keys:
                continue
            residue_tokens.add((chain, resid, resname))
            try:
                center_values.append(
                    (
                        float(line[30:38]),
                        float(line[38:46]),
                        float(line[46:54]),
                    )
                )
            except ValueError:
                continue
    if not residue_tokens:
        return None
    residue_expressions = tuple(
        f"chain {chain} and resid {resid}" for chain, resid, _ in sorted(residue_tokens)
    )
    residue_ids = tuple(
        f"residue:{chain}:{resid}:{resname}" for chain, resid, resname in sorted(residue_tokens)
    )
    center = _mean_center(center_values)
    pocket_id = f"pocket:fpocket:{rank}"
    annotations = {
        "definition_id": "pocket",
        "definition_source": definition("pocket").source,
        "reference_keys": list(definition("pocket").reference_keys),
        "backend": FPOCKET_BACKEND,
        "method": "fpocket_alpha_sphere_geometric_pocket",
        "rank": rank,
        "center_A": list(center) if center is not None else None,
        "lining_residue_count": len(residue_expressions),
        "limitation": "fpocket pockets are geometric cavities, not calibrated ligand-binding-site probabilities.",
    }
    return PocketRecord(
        id=pocket_id,
        name=f"pocket{rank}_fpocket",
        backend=FPOCKET_BACKEND,
        rank=rank,
        center_A=center,
        score=None,
        probability=None,
        residue_ids=residue_ids,
        residue_expressions=residue_expressions,
        annotations=annotations,
    )


def _hydrate_pocket_store(
    store: PocketStore,
    residues_by_key: dict[tuple[str, str], Any],
) -> PocketStore:
    records: list[PocketRecord] = []
    for record in store.records:
        residue_ids: list[str] = []
        expressions: list[str] = []
        for expression in record.residue_expressions:
            chain, resid = _chain_resid_from_expression(expression)
            residue = residues_by_key.get((chain, resid))
            if residue is None:
                continue
            residue_ids.append(object_id_for_residue(residue))
            expressions.append(_residue_expression(residue))
        records.append(
            PocketRecord(
                id=record.id,
                name=record.name,
                backend=record.backend,
                rank=record.rank,
                center_A=record.center_A,
                score=record.score,
                probability=record.probability,
                residue_ids=tuple(dict.fromkeys(residue_ids)),
                residue_expressions=tuple(dict.fromkeys(expressions)),
                annotations={
                    **record.annotations,
                    "lining_residues": list(dict.fromkeys(residue_ids)),
                    "lining_residue_count": len(set(residue_ids)),
                },
            )
        )
    return PocketStore(
        backend=store.backend,
        records=tuple(record for record in records if record.residue_expressions),
        limitations=store.limitations,
    )


def object_ref_for_pocket(record: PocketRecord) -> ObjectRef:
    """Return a public object reference for one backend pocket."""

    return ObjectRef(
        id=record.id,
        type="pocket",
        name=record.name,
        annotations=record.annotations,
    )


def selection_expression_for_pocket(record: PocketRecord) -> str:
    """Return an MDAnalysis selection expression for pocket lining residues."""

    return " or ".join(f"({expression})" for expression in record.residue_expressions)


def _p2rank_command() -> str | None:
    env_command = os.environ.get("MOLINSPECT_P2RANK_COMMAND")
    if env_command:
        return env_command
    return shutil.which("prank") or shutil.which("p2rank")


def _residues_by_chain_resid(universe: Any) -> dict[tuple[str, str], Any]:
    return {
        (chain_for_residue(residue), str(getattr(residue, "resid", "")).strip()): residue
        for residue in universe.residues
    }


def _p2rank_residue_expression(token: str) -> str | None:
    match = _P2RANK_RESIDUE_PATTERN.fullmatch(token.strip())
    if match is None:
        return None
    return f"chain {match.group('chain')} and resid {match.group('resid')}"


def _chain_resid_from_expression(expression: str) -> tuple[str, str]:
    parts = expression.split()
    try:
        return parts[1], parts[4]
    except IndexError:
        return "_", "_"


def _residue_expression(residue: Any) -> str:
    chain = chain_for_residue(residue)
    resid = str(getattr(residue, "resid", "")).strip()
    resname = str(getattr(residue, "resname", "")).strip()
    parts = [f"chain {chain}", f"resid {resid}"]
    if resname:
        parts.append(f"resname {resname}")
    return " and ".join(parts)


def _row_value(row: dict[str, str], key: str) -> str | None:
    for raw_key, raw_value in row.items():
        if raw_key.strip() == key:
            return raw_value.strip()
    return None


def _center_from_row(row: dict[str, str]) -> tuple[float, float, float] | None:
    x = _as_float(_row_value(row, "center_x"))
    y = _as_float(_row_value(row, "center_y"))
    z = _as_float(_row_value(row, "center_z"))
    if x is None or y is None or z is None:
        return None
    return (x, y, z)


def _mean_center(values: list[tuple[float, float, float]]) -> tuple[float, float, float] | None:
    if not values:
        return None
    count = float(len(values))
    return (
        round(sum(value[0] for value in values) / count, 3),
        round(sum(value[1] for value in values) / count, 3),
        round(sum(value[2] for value in values) / count, 3),
    )


def _as_float(value: str | None) -> float | None:
    if value is None or value == "":
        return None
    try:
        return float(value)
    except ValueError:
        return None


def _as_int(value: str | None) -> int | None:
    if value is None or value == "":
        return None
    try:
        return int(value)
    except ValueError:
        return None
