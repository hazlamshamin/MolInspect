"""Temporal analysis backend adapters."""

from __future__ import annotations

import warnings
from collections import Counter
from dataclasses import dataclass
from typing import Any

import numpy as np

from ..bonds import atom_element
from ..definitions import (
    TEMPORAL_HBOND_DONOR_ACCEPTOR_DISTANCE_A,
    TEMPORAL_HBOND_DONOR_HYDROGEN_DISTANCE_A,
    TEMPORAL_HBOND_MIN_ANGLE_DEG,
    definition,
)

MDANALYSIS_HBOND_BACKEND = "MDAnalysis HydrogenBondAnalysis"
DONOR_ACCEPTOR_ELEMENTS = frozenset({"N", "O", "S"})


@dataclass(frozen=True, slots=True)
class HydrogenBondObservation:
    """One frame-level hydrogen bond from MDAnalysis HydrogenBondAnalysis."""

    frame: int
    donor_index: int
    hydrogen_index: int
    acceptor_index: int
    distance_A: float
    angle_deg: float

    @property
    def atom_pair_key(self) -> tuple[int, int]:
        return (self.donor_index, self.acceptor_index)


@dataclass(frozen=True, slots=True)
class HydrogenBondSeries:
    """Hydrogen-bond observations and backend limitations for selected frames."""

    frames: tuple[int, ...]
    observations: tuple[HydrogenBondObservation, ...]
    limitations: tuple[str, ...] = ()

    def observations_by_frame(self) -> dict[int, list[HydrogenBondObservation]]:
        grouped: dict[int, list[HydrogenBondObservation]] = {frame: [] for frame in self.frames}
        for observation in self.observations:
            grouped.setdefault(observation.frame, []).append(observation)
        return grouped

    def counts_by_frame(self) -> list[int]:
        grouped = self.observations_by_frame()
        return [len(grouped.get(frame, ())) for frame in self.frames]

    def best_by_frame(self) -> dict[int, HydrogenBondObservation]:
        grouped = self.observations_by_frame()
        best: dict[int, HydrogenBondObservation] = {}
        for frame, observations in grouped.items():
            if observations:
                best[frame] = sorted(
                    observations,
                    key=lambda item: (-item.angle_deg, item.distance_A, item.donor_index),
                )[0]
        return best

    def pair_occupancies(self) -> list[dict[str, Any]]:
        if not self.frames:
            return []
        counts = Counter(observation.atom_pair_key for observation in self.observations)
        return [
            {
                "donor_index": donor_index,
                "acceptor_index": acceptor_index,
                "frame_count": count,
                "occupancy": round(count / len(self.frames), 3),
            }
            for (donor_index, acceptor_index), count in counts.most_common()
        ]


def hydrogen_bond_series(
    universe: Any,
    group1: Any,
    group2: Any,
    frames: list[int],
    report_missing_requirements: bool = False,
) -> HydrogenBondSeries:
    """Return explicit-hydrogen H-bond observations between two atom groups."""

    if not frames:
        return HydrogenBondSeries(frames=(), observations=())

    group1_indices = _atom_indices(group1)
    group2_indices = _atom_indices(group2)
    combined_indices = tuple(dict.fromkeys((*group1_indices, *group2_indices)))
    donors = _indices_with_elements(universe, combined_indices, DONOR_ACCEPTOR_ELEMENTS)
    acceptors = donors
    hydrogens = _hydrogen_indices(universe, combined_indices)
    missing_limitations = _missing_requirement_limitations(donors, acceptors, hydrogens)
    if missing_limitations:
        return HydrogenBondSeries(
            frames=tuple(frames),
            observations=(),
            limitations=missing_limitations if report_missing_requirements else (),
        )

    try:
        from MDAnalysis.analysis.hydrogenbonds.hbond_analysis import (
            HydrogenBondAnalysis,
        )
    except Exception as exc:
        limitation = f"{MDANALYSIS_HBOND_BACKEND} is unavailable: {exc}"
        return HydrogenBondSeries(
            frames=tuple(frames),
            observations=(),
            limitations=(limitation,) if report_missing_requirements else (),
        )

    try:
        analysis = HydrogenBondAnalysis(
            universe=universe,
            donors_sel=_index_selection(donors),
            hydrogens_sel=_index_selection(hydrogens),
            acceptors_sel=_index_selection(acceptors),
            between=[_index_selection(group1_indices), _index_selection(group2_indices)],
            d_h_cutoff=TEMPORAL_HBOND_DONOR_HYDROGEN_DISTANCE_A,
            d_a_cutoff=TEMPORAL_HBOND_DONOR_ACCEPTOR_DISTANCE_A,
            d_h_a_angle_cutoff=TEMPORAL_HBOND_MIN_ANGLE_DEG,
            update_selections=True,
        )
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            analysis.run(frames=frames)
    except Exception as exc:
        limitation = f"{MDANALYSIS_HBOND_BACKEND} failed: {exc}"
        return HydrogenBondSeries(
            frames=tuple(frames),
            observations=(),
            limitations=(limitation,) if report_missing_requirements else (),
        )

    return HydrogenBondSeries(
        frames=tuple(frames),
        observations=tuple(_observations_from_hba_results(analysis.results.hbonds)),
        limitations=(),
    )


def _observations_from_hba_results(results: Any) -> list[HydrogenBondObservation]:
    rows = np.asarray(results)
    if rows.size == 0:
        return []
    observations: list[HydrogenBondObservation] = []
    for row in rows:
        observations.append(
            HydrogenBondObservation(
                frame=int(row[0]),
                donor_index=int(row[1]),
                hydrogen_index=int(row[2]),
                acceptor_index=int(row[3]),
                distance_A=round(float(row[4]), 3),
                angle_deg=round(float(row[5]), 1),
            )
        )
    return observations


def _missing_requirement_limitations(
    donors: tuple[int, ...],
    acceptors: tuple[int, ...],
    hydrogens: tuple[int, ...],
) -> tuple[str, ...]:
    limitations: list[str] = []
    if not donors or not acceptors:
        limitations.append(
            "HydrogenBondAnalysis requires selected donor/acceptor atoms with N, O, or S elements."
        )
    if not hydrogens:
        limitations.append(definition("hydrogen_bond_occupancy").limitations[0])
    return tuple(limitations)


def _atom_indices(atomgroup: Any) -> tuple[int, ...]:
    return tuple(int(atom.ix) for atom in atomgroup)


def _indices_with_elements(
    universe: Any,
    indices: tuple[int, ...],
    elements: frozenset[str],
) -> tuple[int, ...]:
    return tuple(index for index in indices if atom_element(universe.atoms[index]) in elements)


def _hydrogen_indices(universe: Any, indices: tuple[int, ...]) -> tuple[int, ...]:
    hydrogen_indices: list[int] = []
    for index in indices:
        atom = universe.atoms[index]
        element = atom_element(atom)
        name = str(getattr(atom, "name", "")).strip().upper()
        if element == "H" or name.startswith("H"):
            hydrogen_indices.append(index)
    return tuple(hydrogen_indices)


def _index_selection(indices: tuple[int, ...]) -> str:
    return "index " + " ".join(str(index) for index in indices)
