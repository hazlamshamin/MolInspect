"""Geometry and timeline metric helpers."""

from __future__ import annotations

from collections.abc import Iterable
from dataclasses import dataclass
from typing import Any

import numpy as np

from .errors import MetricError
from .notation import FRAME_RANGE_NOTATION_HELP
from .selections import normalize_frame
from .definitions import CONTACT_CUTOFF_A

MAX_SAMPLED_VALUES = 50


@dataclass(frozen=True, slots=True)
class AtomPairDistance:
    """Closest heavy-atom pair between two atom groups."""

    distance_A: float
    atom_a: Any
    atom_b: Any


@dataclass(frozen=True, slots=True)
class FitTransform:
    """Rigid transform that superposes mobile coordinates onto reference coordinates."""

    reference_center: np.ndarray
    mobile_center: np.ndarray
    rotation: np.ndarray


def center_of_geometry_A(atomgroup: Any) -> list[float]:
    """Return center of geometry in angstroms, rounded for compact output."""

    if len(atomgroup) == 0:
        raise MetricError("cannot compute center of geometry for an empty selection")
    center = np.asarray(atomgroup.positions, dtype=float).mean(axis=0)
    return _round_vector(center)


def min_distance_A(atomgroup_a: Any, atomgroup_b: Any, heavy_only: bool = True) -> float | None:
    """Return minimum pairwise distance in angstroms."""

    group_a = heavy_atomgroup(atomgroup_a) if heavy_only else atomgroup_a
    group_b = heavy_atomgroup(atomgroup_b) if heavy_only else atomgroup_b
    if len(group_a) == 0 or len(group_b) == 0:
        return None

    distances = distance_matrix_A(group_a, group_b)
    return round(float(np.min(distances)), 3)


def closest_heavy_atom_pair(atomgroup_a: Any, atomgroup_b: Any) -> AtomPairDistance | None:
    """Return the closest heavy-atom pair between two atom groups."""

    group_a = heavy_atomgroup(atomgroup_a)
    group_b = heavy_atomgroup(atomgroup_b)
    if len(group_a) == 0 or len(group_b) == 0:
        return None

    atoms_a = list(group_a)
    atoms_b = list(group_b)
    distances = distance_matrix_A(group_a, group_b)
    flat_index = int(np.argmin(distances))
    index_a, index_b = np.unravel_index(flat_index, distances.shape)
    return AtomPairDistance(
        distance_A=round(float(distances[index_a, index_b]), 3),
        atom_a=atoms_a[index_a],
        atom_b=atoms_b[index_b],
    )


def distance_between_centers_A(atomgroup_a: Any, atomgroup_b: Any) -> float | None:
    """Return distance between centers of geometry in angstroms."""

    if len(atomgroup_a) == 0 or len(atomgroup_b) == 0:
        return None
    center_a = np.asarray(atomgroup_a.positions, dtype=float).mean(axis=0)
    center_b = np.asarray(atomgroup_b.positions, dtype=float).mean(axis=0)
    return round(float(np.linalg.norm(center_a - center_b)), 3)


def aligned_rmsd_A(reference_positions: Any, mobile_positions: Any) -> float:
    """Return Kabsch-aligned RMSD between two coordinate sets."""

    reference = _as_position_array(reference_positions)
    mobile = _as_position_array(mobile_positions)
    _validate_position_pair(reference, mobile)
    try:
        from MDAnalysis.analysis.rms import rmsd as mdanalysis_rmsd

        return round(
            float(mdanalysis_rmsd(mobile, reference, center=True, superposition=True)),
            3,
        )
    except Exception:
        pass
    aligned = align_positions_to_reference(reference, mobile)
    deltas = aligned - reference
    return round(float(np.sqrt(np.mean(np.sum(deltas * deltas, axis=1)))), 3)


def align_positions_to_reference(reference_positions: Any, mobile_positions: Any) -> np.ndarray:
    """Return mobile coordinates optimally superposed on reference coordinates."""

    transform = fit_transform_to_reference(reference_positions, mobile_positions)
    return apply_fit_transform(mobile_positions, transform)


def fit_transform_to_reference(reference_positions: Any, mobile_positions: Any) -> FitTransform:
    """Return the Kabsch transform that superposes mobile coordinates on reference."""

    reference = _as_position_array(reference_positions)
    mobile = _as_position_array(mobile_positions)
    _validate_position_pair(reference, mobile)

    reference_center = reference.mean(axis=0)
    mobile_center = mobile.mean(axis=0)
    reference_centered = reference - reference_center
    mobile_centered = mobile - mobile_center
    covariance = mobile_centered.T @ reference_centered
    left, _, right_t = np.linalg.svd(covariance)
    correction = np.eye(3)
    correction[2, 2] = np.sign(np.linalg.det(left @ right_t)) or 1.0
    rotation = left @ correction @ right_t
    return FitTransform(
        reference_center=reference_center,
        mobile_center=mobile_center,
        rotation=rotation,
    )


def apply_fit_transform(positions: Any, transform: FitTransform) -> np.ndarray:
    """Apply a Kabsch transform from another atom set to arbitrary coordinates."""

    coordinates = _as_position_array(positions)
    return (coordinates - transform.mobile_center) @ transform.rotation + transform.reference_center


def distance_matrix_A(atomgroup_a: Any, atomgroup_b: Any) -> np.ndarray:
    """Return pairwise distances in angstroms, using MDAnalysis PBC support when available."""

    positions_a = np.asarray(atomgroup_a.positions, dtype=float)
    positions_b = np.asarray(atomgroup_b.positions, dtype=float)
    try:
        from MDAnalysis.lib.distances import distance_array

        return np.asarray(
            distance_array(positions_a, positions_b, box=_current_box(atomgroup_a)),
            dtype=float,
        )
    except Exception:
        deltas = positions_a[:, None, :] - positions_b[None, :, :]
        return np.sqrt(np.sum(deltas * deltas, axis=2))


def rmsf_per_atom_A(position_frames: list[Any], align: bool = True) -> list[float]:
    """Return per-atom RMSF across frames, optionally aligned to the first frame."""

    if not position_frames:
        raise MetricError("cannot compute RMSF without selected frames")
    reference = _as_position_array(position_frames[0])
    if align:
        positions = [align_positions_to_reference(reference, frame) for frame in position_frames]
    else:
        positions = [_as_position_array(frame) for frame in position_frames]
    stack = np.stack(positions)
    try:
        import MDAnalysis as mda
        from MDAnalysis.analysis.rms import RMSF

        universe = mda.Universe.empty(stack.shape[1], trajectory=True)
        universe.load_new(stack, order="fac")
        results = RMSF(universe.atoms).run()
        return [round(float(value), 3) for value in results.results.rmsf.tolist()]
    except Exception:
        pass
    mean_positions = stack.mean(axis=0)
    deltas = stack - mean_positions
    values = np.sqrt(np.mean(np.sum(deltas * deltas, axis=2), axis=0))
    return [round(float(value), 3) for value in values.tolist()]


def heavy_atomgroup(atomgroup: Any) -> Any:
    """Return heavy atoms, falling back to the original group if metadata is sparse."""

    heavy_indices = [atom.ix for atom in atomgroup if not _is_hydrogen(atom)]
    if not heavy_indices:
        return atomgroup
    return atomgroup.universe.atoms[heavy_indices]


def relation_type(distance_A: float | None, contact_cutoff_A: float = CONTACT_CUTOFF_A) -> str:
    if distance_A is not None and distance_A <= contact_cutoff_A:
        return "nonbonded_contact"
    return "near"


def frame_indices(n_frames: int, frames: str | tuple[int, int] = "all", stride: int | None = None) -> list[int]:
    """Resolve a public frame range into explicit frame indices."""

    if stride is not None and stride < 1:
        raise MetricError("stride must be >= 1")
    step = stride or 1

    if frames == "all":
        return list(range(0, n_frames, step))
    if isinstance(frames, tuple) and len(frames) == 2:
        start = normalize_frame(frames[0], n_frames)
        end = normalize_frame(frames[1], n_frames)
        if end < start:
            raise MetricError(f"frame range end must be >= start. {FRAME_RANGE_NOTATION_HELP}")
        return list(range(start, end + 1, step))
    raise MetricError(FRAME_RANGE_NOTATION_HELP)


def summarize_distances(
    distances_A: list[float | None],
    contact_cutoff_A: float = CONTACT_CUTOFF_A,
) -> dict[str, Any]:
    numeric = np.array([value for value in distances_A if value is not None], dtype=float)
    if numeric.size == 0:
        return {
            "min_A": None,
            "median_A": None,
            "max_A": None,
            "mean_A": None,
            "contact_cutoff_A": contact_cutoff_A,
            "contact_occupancy": 0.0,
        }

    contact_occupancy = float(np.mean(numeric <= contact_cutoff_A))
    return {
        "min_A": round(float(np.min(numeric)), 3),
        "median_A": round(float(np.median(numeric)), 3),
        "max_A": round(float(np.max(numeric)), 3),
        "mean_A": round(float(np.mean(numeric)), 3),
        "contact_cutoff_A": contact_cutoff_A,
        "contact_occupancy": round(contact_occupancy, 3),
    }


def contact_events(
    frames: list[int],
    distances_A: list[float | None],
    contact_cutoff_A: float = CONTACT_CUTOFF_A,
) -> list[dict[str, Any]]:
    """Return simple contact formation/break events."""

    events: list[dict[str, Any]] = []
    previous: bool | None = None
    for frame, distance in zip(frames, distances_A, strict=True):
        if distance is None:
            current = False
        else:
            current = distance <= contact_cutoff_A
        if previous is not None and current != previous:
            events.append(
                {
                    "type": "contact_forms" if current else "contact_breaks",
                    "frame": frame,
                    "distance_A": distance,
                    "cutoff_A": contact_cutoff_A,
                }
            )
        previous = current
    return events


def sampled_values(
    frames: list[int],
    values: Iterable[float | bool | None],
    value_key: str,
    max_samples: int = MAX_SAMPLED_VALUES,
) -> list[dict[str, Any]]:
    """Return compact frame-sampled values."""

    value_list = list(values)
    if len(frames) <= max_samples:
        chosen = list(range(len(frames)))
    else:
        chosen = sorted(set(np.linspace(0, len(frames) - 1, max_samples, dtype=int).tolist()))
    return [{"frame": frames[index], value_key: value_list[index]} for index in chosen]


def _is_hydrogen(atom: Any) -> bool:
    element = str(getattr(atom, "element", "")).strip().upper()
    if element == "H":
        return True
    name = str(getattr(atom, "name", "")).strip().upper()
    return name.startswith("H")


def _round_vector(vector: np.ndarray) -> list[float]:
    return [round(float(value), 3) for value in vector.tolist()]


def _as_position_array(positions: Any) -> np.ndarray:
    array = np.asarray(positions, dtype=float)
    if array.ndim != 2 or array.shape[1] != 3:
        raise MetricError("coordinate arrays must have shape (n_atoms, 3)")
    return array


def _validate_position_pair(reference: np.ndarray, mobile: np.ndarray) -> None:
    if reference.shape != mobile.shape:
        raise MetricError(
            "RMSD/RMSF selections must resolve to the same atom count in every frame."
        )
    if len(reference) == 0:
        raise MetricError("cannot align an empty selection")


def _current_box(atomgroup: Any) -> Any | None:
    dimensions = getattr(getattr(atomgroup, "universe", None), "dimensions", None)
    if dimensions is None:
        return None
    array = np.asarray(dimensions, dtype=float)
    if array.shape != (6,) or not np.all(np.isfinite(array)) or np.any(array[:3] <= 0):
        return None
    return array
