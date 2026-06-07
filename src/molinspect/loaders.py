"""Structure and trajectory loading."""

from __future__ import annotations

import hashlib
import tempfile
from pathlib import Path
from typing import Any

from .backends import available_backend_annotation_ids
from .errors import LoadError
from .objects import chain_identifier_for_atom
from .schemas import LoadResult
from .session import InspectionSession
from .world import InspectionWorld


def load(
    structure: str | Path | None = None,
    topology: str | Path | None = None,
    trajectory: str | Path | None = None,
    name: str | None = None,
) -> InspectionSession:
    """Load a static structure or trajectory into an inspection session."""

    structure_path = _optional_path(structure)
    topology_path = _optional_path(topology)
    trajectory_path = _optional_path(trajectory)

    if structure_path is not None and (topology_path is not None or trajectory_path is not None):
        raise LoadError("Use either structure=... or topology=/trajectory=..., not both.")
    if structure_path is None and topology_path is None:
        raise LoadError("Provide structure=... for static input or topology=... for trajectory input.")
    if trajectory_path is not None and topology_path is None:
        raise LoadError("trajectory=... requires topology=...")

    source_files = tuple(
        path for path in (structure_path, topology_path, trajectory_path) if path is not None
    )
    for path in source_files:
        if not path.exists():
            raise LoadError(f"Input file does not exist: {path}")
        if not path.is_file():
            raise LoadError(f"Input path is not a file: {path}")

    universe = _load_universe(structure_path, topology_path, trajectory_path)
    target_id = _target_id(source_files)
    load_name = name or _default_name(source_files)
    load_result = _load_result(universe, target_id, load_name, source_files)
    world = InspectionWorld(universe=universe, source_files=source_files)

    return InspectionSession(
        universe=universe,
        world=world,
        target_id=target_id,
        name=load_name,
        source_files=source_files,
        load_result=load_result,
    )


def _load_universe(
    structure: Path | None,
    topology: Path | None,
    trajectory: Path | None,
) -> Any:
    try:
        import MDAnalysis as mda
    except ImportError as exc:
        raise LoadError("MDAnalysis is required for the default loader. Install with `uv sync`.") from exc

    try:
        if structure is not None:
            if structure.suffix.lower() in {".cif", ".mmcif"}:
                return _load_mmcif_with_gemmi(mda, structure)
            return mda.Universe(str(structure))
        if trajectory is not None:
            return mda.Universe(str(topology), str(trajectory))
        return mda.Universe(str(topology))
    except Exception as exc:
        raise LoadError(f"MDAnalysis could not load the supplied input: {exc}") from exc


def _load_mmcif_with_gemmi(mda: Any, structure: Path) -> Any:
    """Convert mmCIF to an in-memory MDAnalysis universe through optional Gemmi."""

    try:
        import gemmi  # type: ignore[import-not-found]
    except ImportError as exc:
        raise LoadError(
            "mmCIF input requires optional Gemmi support. Install with `uv sync --extra static`."
        ) from exc

    try:
        gemmi_structure = gemmi.read_structure(str(structure))
        with tempfile.NamedTemporaryFile(suffix=".pdb") as converted:
            gemmi_structure.write_pdb(converted.name)
            return mda.Universe(converted.name, in_memory=True)
    except Exception as exc:
        raise LoadError(f"Gemmi could not convert mmCIF input to PDB for loading: {exc}") from exc


def _load_result(
    universe: Any,
    target_id: str,
    name: str | None,
    source_files: tuple[Path, ...],
) -> LoadResult:
    n_frames = len(universe.trajectory)
    return LoadResult(
        target_id=target_id,
        name=name,
        n_atoms=len(universe.atoms),
        n_residues=len(universe.residues),
        n_chains=_count_chains(universe),
        n_segments=len(universe.segments),
        n_frames=n_frames,
        mode="static" if n_frames == 1 else "trajectory",
        source_files=[str(path) for path in source_files],
        available_annotations=_available_annotations(source_files),
    )


def _count_chains(universe: Any) -> int:
    chains = {chain_identifier_for_atom(atom) for atom in universe.atoms}
    return len(chains)


def _available_annotations(source_files: tuple[Path, ...]) -> list[str]:
    return available_backend_annotation_ids(source_files)


def _target_id(source_files: tuple[Path, ...]) -> str:
    digest = hashlib.sha256()
    digest.update(b"molinspect-target-v1")
    for path in source_files:
        resolved = path.resolve()
        digest.update(str(resolved).encode())
        digest.update(b"\0")
        with resolved.open("rb") as handle:
            for chunk in iter(lambda: handle.read(1024 * 1024), b""):
                digest.update(chunk)
        digest.update(b"\0")
    return f"sha256:{digest.hexdigest()}"


def _optional_path(value: str | Path | None) -> Path | None:
    if value is None:
        return None
    return Path(value).expanduser()


def _default_name(source_files: tuple[Path, ...]) -> str | None:
    if not source_files:
        return None
    return source_files[0].stem
