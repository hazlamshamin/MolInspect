"""Temporal inspection primitives for trajectories and frame series."""

from __future__ import annotations

from typing import Any

import numpy as np

from ..backends.temporal import (
    MDANALYSIS_HBOND_BACKEND,
    HydrogenBondObservation,
    HydrogenBondSeries,
    hydrogen_bond_series,
)
from ..definitions import (
    CONTACT_CUTOFF_A,
    LIGAND_STABILITY_STABLE_MAX_RMSD_A,
    LIGAND_STABILITY_STABLE_MIN_CONTACT_OCCUPANCY,
    LIGAND_STABILITY_UNSTABLE_MAX_CONTACT_OCCUPANCY,
    LIGAND_STABILITY_UNSTABLE_MIN_RMSD_A,
    STATE_MIN_RMSD_RANGE_A,
    TEMPORAL_HBOND_DONOR_ACCEPTOR_DISTANCE_A,
    TEMPORAL_HBOND_DONOR_HYDROGEN_DISTANCE_A,
    TEMPORAL_HBOND_MIN_ANGLE_DEG,
    definition,
)
from ..errors import MetricError
from ..metrics import (
    apply_fit_transform,
    align_positions_to_reference,
    aligned_rmsd_A,
    contact_events,
    distance_between_centers_A,
    fit_transform_to_reference,
    frame_indices,
    rmsf_per_atom_A,
    sampled_values,
    summarize_distances,
)
from ..notation import TIMELINE_METRIC_NOTATION_HELP
from ..objects import object_id_for_residue
from ..relations import atom_label, relation_priority_for_type
from ..schemas import EvidenceItem, Relation, TimelineResult
from ..world import InspectionWorld
from .events import (
    centroid_distance_events,
    movement_events,
    numeric_summary,
    relation_change_events,
    representative_displacement_frames,
    representative_distance_frames,
    representative_numeric_frames,
    representative_relation_frames,
    selection_spread_events,
    selection_spread_representative_frames,
)
from .summaries import relation_counts
from .timeline_specs import SUPPORTED_TIMELINE_METRICS


def timeline(
    world: InspectionWorld,
    metric: str,
    selection: str | None = None,
    selection1: str | None = None,
    selection2: str | None = None,
    frames: str | tuple[int, int] = "all",
    stride: int | None = None,
) -> TimelineResult:
    """Summarize structural metrics, relations, or mobility over frames."""

    universe = world.universe
    if not isinstance(metric, str):
        raise MetricError(TIMELINE_METRIC_NOTATION_HELP)
    metric_name = metric.strip().lower()
    if metric_name not in SUPPORTED_TIMELINE_METRICS:
        raise MetricError(TIMELINE_METRIC_NOTATION_HELP)

    selected_frames = frame_indices(len(universe.trajectory), frames=frames, stride=stride)
    if not selected_frames:
        raise MetricError("no frames selected for timeline")

    if metric_name == "rmsd":
        if selection is None:
            raise MetricError("timeline metric='rmsd' requires selection.")
        return _rmsd_timeline(world, selection, selected_frames)
    if metric_name == "states":
        if selection is None:
            raise MetricError("timeline metric='states' requires selection.")
        return _states_timeline(world, selection, selected_frames)
    if metric_name in {"rmsf", "mobility"}:
        if selection is None:
            raise MetricError("timeline metric='rmsf' or metric='mobility' requires selection.")
        return _mobility_timeline(world, metric_name, selection, selected_frames)
    if metric_name == "displacement":
        if selection is None:
            raise MetricError("timeline metric='displacement' requires selection.")
        return _displacement_timeline(world, selection, selected_frames)
    if metric_name == "selection_spread":
        if selection is None:
            raise MetricError("timeline metric='selection_spread' requires selection.")
        return _selection_spread_timeline(world, selection, selected_frames)

    if selection1 is None or selection2 is None:
        raise MetricError(
            "timeline requires selection1 and selection2 for metric='distance', "
            "metric='contact', metric='relation', metric='hydrogen_bonds', "
            "metric='interaction_persistence', metric='ligand_stability', "
            "and metric='centroid_distance'."
        )

    resolved1 = world.resolve_selection(selection1, frame=selected_frames[0])
    resolved2 = world.resolve_selection(selection2, frame=selected_frames[0])
    if metric_name == "centroid_distance":
        return _centroid_distance_timeline(
            world,
            selection1,
            selection2,
            selected_frames,
            resolved1,
            resolved2,
        )
    if metric_name == "hydrogen_bonds":
        return _hydrogen_bond_timeline(
            world,
            selection1,
            selection2,
            selected_frames,
            resolved1,
            resolved2,
        )
    if metric_name == "interaction_persistence":
        return _interaction_persistence_timeline(
            world,
            selection1,
            selection2,
            selected_frames,
            resolved1,
            resolved2,
        )
    if metric_name == "ligand_stability":
        return _ligand_stability_timeline(
            world,
            selection1,
            selection2,
            selected_frames,
            resolved1,
            resolved2,
        )

    return _relation_timeline_core(
        world,
        selection1,
        selection2,
        selected_frames,
        resolved1,
        resolved2,
        metric_name=metric_name,
    )


def _centroid_distance_timeline(
    world: InspectionWorld,
    selection1: str,
    selection2: str,
    selected_frames: list[int],
    resolved1: Any,
    resolved2: Any,
) -> TimelineResult:
    universe = world.universe
    distances: list[float | None] = []
    for frame in selected_frames:
        universe.trajectory[frame]
        distances.append(distance_between_centers_A(resolved1.atomgroup, resolved2.atomgroup))

    summary = numeric_summary(
        [distance for distance in distances if distance is not None],
        "centroid_distance_A",
    )
    events = centroid_distance_events(selected_frames, distances)
    representatives = representative_numeric_frames(
        selected_frames,
        [float(distance) for distance in distances if distance is not None],
        "centroid_distance_A",
    )
    if "lowest" in representatives:
        representatives["closest_centroids"] = representatives.pop("lowest")
    if "highest" in representatives:
        representatives["farthest_centroids"] = representatives.pop("highest")
    summary.update({"center": "center_of_geometry", "representative_frames": representatives})

    return _timeline_result(
        metric="centroid_distance",
        selection1=selection1,
        selection2=selection2,
        frames_analyzed=len(selected_frames),
        summary=summary,
        events=events,
        sampled_values=sampled_values(selected_frames, distances, "centroid_distance_A"),
        limitations=list(
            dict.fromkeys(
                [
                    definition("centroid_distance").description,
                    *definition("centroid_distance").limitations,
                    *resolved1.limitations,
                    *resolved2.limitations,
                ]
            )
        ),
    )


def _relation_timeline_core(
    world: InspectionWorld,
    selection1: str,
    selection2: str,
    selected_frames: list[int],
    resolved1: Any,
    resolved2: Any,
    metric_name: str,
) -> TimelineResult:
    universe = world.universe
    source_id = _timeline_source_id(resolved1.resolved_objects, selection1)
    target_id = _timeline_source_id(resolved2.resolved_objects, selection2)
    hbond_series = hydrogen_bond_series(
        universe,
        resolved1.atomgroup,
        resolved2.atomgroup,
        selected_frames,
    )
    hbonds_by_frame = hbond_series.best_by_frame()
    relations: list[Relation | None] = []
    distances: list[float | None] = []
    for frame in selected_frames:
        universe.trajectory[frame]
        hbond = hbonds_by_frame.get(frame)
        relation: Relation | None
        if hbond is not None:
            relation = _relation_from_hbond_observation(
                universe,
                hbond,
                resolved1.atomgroup,
                resolved2.atomgroup,
                source_id,
                target_id,
            )
        else:
            relation = world.relation_between(
                resolved1.atomgroup,
                resolved2.atomgroup,
                source_id=source_id,
                target_id=target_id,
            )
        relations.append(relation)
        distances.append(relation.min_distance_A if relation is not None else None)

    summary = summarize_distances(distances)
    relation_type_counts = relation_counts([relation for relation in relations if relation is not None])
    if relation_type_counts:
        summary["relation_type_counts"] = relation_type_counts
        summary["relation_occupancy"] = {
            relation_type: round(count / len(selected_frames), 3)
            for relation_type, count in relation_type_counts.items()
        }
        summary["dominant_relation_type"] = max(
            relation_type_counts.items(),
            key=lambda item: (item[1], -_relation_priority(item[0])),
        )[0]
        summary["representative_frames"] = representative_relation_frames(
            selected_frames,
            relations,
            distances,
        )
    events = contact_events(selected_frames, distances)
    if metric_name == "relation":
        events.extend(relation_change_events(selected_frames, relations))
    events.extend(movement_events(selected_frames, distances, source_id, target_id))
    contact_values = [
        None if distance is None else distance <= CONTACT_CUTOFF_A for distance in distances
    ]
    if metric_name == "distance":
        samples = sampled_values(selected_frames, distances, "value_A")
    elif metric_name == "contact":
        samples = sampled_values(selected_frames, contact_values, "contact")
    else:
        samples = _sampled_relations(selected_frames, relations)
    summary.setdefault(
        "representative_frames",
        representative_distance_frames(selected_frames, distances, events),
    )

    limitations = [
        "Distances are minimum heavy-atom distances using MDAnalysis distance arrays when available.",
        f"Contacts use a <= {CONTACT_CUTOFF_A:g} A cutoff.",
    ]
    if hbond_series.observations:
        limitations.append(
            "Hydrogen-bond relation frames use MDAnalysis HydrogenBondAnalysis; other relation "
            "types use closest heavy-atom heuristic rules."
        )
    else:
        limitations.append(
            "Relation types use closest heavy-atom heuristic rules when no explicit-hydrogen "
            "hydrogen-bond backend evidence is available."
        )
    limitations.extend(resolved1.limitations)
    limitations.extend(resolved2.limitations)

    return _timeline_result(
        metric=metric_name,
        selection1=selection1,
        selection2=selection2,
        frames_analyzed=len(selected_frames),
        summary=summary,
        events=events,
        sampled_values=samples,
        limitations=list(dict.fromkeys(limitations)),
    )


def _hydrogen_bond_timeline(
    world: InspectionWorld,
    selection1: str,
    selection2: str,
    selected_frames: list[int],
    resolved1: Any,
    resolved2: Any,
) -> TimelineResult:
    universe = world.universe
    source_id = _timeline_source_id(resolved1.resolved_objects, selection1)
    target_id = _timeline_source_id(resolved2.resolved_objects, selection2)
    series = hydrogen_bond_series(
        universe,
        resolved1.atomgroup,
        resolved2.atomgroup,
        selected_frames,
        report_missing_requirements=True,
    )
    counts = series.counts_by_frame()
    present = [count > 0 for count in counts]
    summary = numeric_summary(counts, "hydrogen_bond_count")
    summary.update(
        {
            "backend": MDANALYSIS_HBOND_BACKEND,
            "hydrogen_bond_occupancy": round(sum(present) / len(selected_frames), 3),
            "total_hydrogen_bond_observations": len(series.observations),
            "donor_acceptor_distance_cutoff_A": TEMPORAL_HBOND_DONOR_ACCEPTOR_DISTANCE_A,
            "donor_hydrogen_distance_cutoff_A": TEMPORAL_HBOND_DONOR_HYDROGEN_DISTANCE_A,
            "donor_hydrogen_acceptor_angle_min_deg": TEMPORAL_HBOND_MIN_ANGLE_DEG,
            "representative_frames": _hydrogen_bond_representative_frames(selected_frames, counts),
            "top_hydrogen_bond_pairs": _top_hydrogen_bond_pairs(universe, series),
        }
    )
    limitations = [
        definition("hydrogen_bond_occupancy").description,
        *definition("hydrogen_bond_occupancy").limitations,
        *series.limitations,
        *resolved1.limitations,
        *resolved2.limitations,
    ]
    return _timeline_result(
        metric="hydrogen_bonds",
        selection1=selection1,
        selection2=selection2,
        frames_analyzed=len(selected_frames),
        summary=summary,
        events=_hydrogen_bond_events(selected_frames, counts, source_id, target_id),
        sampled_values=_sampled_hydrogen_bonds(universe, selected_frames, series),
        limitations=list(dict.fromkeys(limitations)),
    )


def _interaction_persistence_timeline(
    world: InspectionWorld,
    selection1: str,
    selection2: str,
    selected_frames: list[int],
    resolved1: Any,
    resolved2: Any,
) -> TimelineResult:
    relation_result = _relation_timeline_core(
        world,
        selection1,
        selection2,
        selected_frames,
        resolved1,
        resolved2,
        metric_name="relation",
    )
    hbond_result = _hydrogen_bond_timeline(
        world,
        selection1,
        selection2,
        selected_frames,
        resolved1,
        resolved2,
    )
    hbond_summary = hbond_result.summary
    relation_summary = relation_result.summary
    contact_occupancy = float(relation_summary.get("contact_occupancy") or 0.0)
    hbond_occupancy = float(hbond_summary.get("hydrogen_bond_occupancy") or 0.0)
    summary = {
        "contact_occupancy": round(contact_occupancy, 3),
        "hydrogen_bond_occupancy": round(hbond_occupancy, 3),
        "dominant_relation_type": relation_summary.get("dominant_relation_type"),
        "relation_occupancy": relation_summary.get("relation_occupancy", {}),
        "min_distance_A": relation_summary.get("min_A"),
        "median_distance_A": relation_summary.get("median_A"),
        "max_distance_A": relation_summary.get("max_A"),
        "contact_cutoff_A": CONTACT_CUTOFF_A,
        "representative_frames": relation_summary.get("representative_frames", {}),
        "top_hydrogen_bond_pairs": hbond_summary.get("top_hydrogen_bond_pairs", []),
    }
    events = [
        *_key_persistence_events(relation_result.events),
        *_key_persistence_events(hbond_result.events),
    ]
    samples = _interaction_persistence_samples(
        relation_result.sampled_values,
        hbond_result.sampled_values,
    )
    limitations = [
        definition("interaction_persistence").description,
        *definition("interaction_persistence").limitations,
        *relation_result.limitations,
        *hbond_result.limitations,
    ]
    return _timeline_result(
        metric="interaction_persistence",
        selection1=selection1,
        selection2=selection2,
        frames_analyzed=len(selected_frames),
        summary=summary,
        events=events,
        sampled_values=samples,
        limitations=list(dict.fromkeys(limitations)),
    )


def _ligand_stability_timeline(
    world: InspectionWorld,
    selection1: str,
    selection2: str,
    selected_frames: list[int],
    resolved1: Any,
    resolved2: Any,
) -> TimelineResult:
    universe = world.universe
    ligand_id = _timeline_source_id(resolved1.resolved_objects, selection1)
    site_id = _timeline_source_id(resolved2.resolved_objects, selection2)
    universe.trajectory[selected_frames[0]]
    ligand_reference = np.asarray(resolved1.atomgroup.positions, dtype=float).copy()
    site_reference = np.asarray(resolved2.atomgroup.positions, dtype=float).copy()

    ligand_rmsd_values: list[float] = []
    center_distances: list[float | None] = []
    for frame in selected_frames:
        universe.trajectory[frame]
        transform = fit_transform_to_reference(site_reference, resolved2.atomgroup.positions)
        aligned_ligand = apply_fit_transform(resolved1.atomgroup.positions, transform)
        ligand_rmsd_values.append(_rmsd_between_position_arrays(ligand_reference, aligned_ligand))
        center_distances.append(distance_between_centers_A(resolved1.atomgroup, resolved2.atomgroup))

    relation_result = _relation_timeline_core(
        world,
        selection1,
        selection2,
        selected_frames,
        resolved1,
        resolved2,
        metric_name="relation",
    )
    hbond_result = _hydrogen_bond_timeline(
        world,
        selection1,
        selection2,
        selected_frames,
        resolved1,
        resolved2,
    )
    contact_occupancy = float(relation_result.summary.get("contact_occupancy") or 0.0)
    hbond_occupancy = float(hbond_result.summary.get("hydrogen_bond_occupancy") or 0.0)
    max_rmsd = max(ligand_rmsd_values) if ligand_rmsd_values else 0.0
    summary = numeric_summary(ligand_rmsd_values, "site_aligned_ligand_rmsd_A")
    center_summary = numeric_summary(
        [distance for distance in center_distances if distance is not None],
        "ligand_site_centroid_distance_A",
    )
    summary.update(center_summary)
    summary.update(
        {
            "ligand_selection": selection1,
            "site_selection": selection2,
            "ligand_object": ligand_id,
            "site_object": site_id,
            "alignment": "binding_site_kabsch",
            "contact_occupancy": round(contact_occupancy, 3),
            "hydrogen_bond_occupancy": round(hbond_occupancy, 3),
            "dominant_relation_type": relation_result.summary.get("dominant_relation_type"),
            "stability_class": _ligand_stability_class(contact_occupancy, max_rmsd),
            "representative_frames": _ligand_stability_representative_frames(
                selected_frames,
                ligand_rmsd_values,
                center_distances,
            ),
        }
    )
    limitations = [
        definition("ligand_stability").description,
        *definition("ligand_stability").limitations,
        *relation_result.limitations,
        *hbond_result.limitations,
    ]
    return _timeline_result(
        metric="ligand_stability",
        selection1=selection1,
        selection2=selection2,
        frames_analyzed=len(selected_frames),
        summary=summary,
        events=_ligand_stability_events(
            selected_frames,
            ligand_rmsd_values,
            relation_result.events,
            ligand_id,
            site_id,
        ),
        sampled_values=_ligand_stability_samples(
            selected_frames,
            ligand_rmsd_values,
            center_distances,
            relation_result.sampled_values,
            hbond_result.sampled_values,
        ),
        limitations=list(dict.fromkeys(limitations)),
    )


def _states_timeline(
    world: InspectionWorld,
    selection: str,
    selected_frames: list[int],
) -> TimelineResult:
    universe = world.universe
    resolved = world.resolve_selection(selection, frame=selected_frames[0])
    universe.trajectory[selected_frames[0]]
    reference_positions = np.asarray(resolved.atomgroup.positions, dtype=float).copy()
    rmsd_values: list[float] = []
    spread_values: list[float] = []
    for frame in selected_frames:
        universe.trajectory[frame]
        rmsd_values.append(aligned_rmsd_A(reference_positions, resolved.atomgroup.positions))
        spread_values.append(_radius_of_gyration_A(resolved.atomgroup.positions))

    assignments = _state_assignments(rmsd_values)
    states = _state_summaries(selected_frames, rmsd_values, spread_values, assignments)
    summary = {
        "state_count": len(states),
        "state_method": "deterministic_rmsd_split",
        "features": ["aligned_rmsd_A", "selection_spread_A"],
        "reference_frame": selected_frames[0],
        "rmsd_range_A": round(max(rmsd_values) - min(rmsd_values), 3) if rmsd_values else 0.0,
        "states": states,
        "representative_frames": {
            state["state_id"]: {
                "frame": state["representative_frame"],
                "mean_rmsd_A": state["mean_rmsd_A"],
            }
            for state in states
        },
    }
    return _timeline_result(
        metric="states",
        selection=selection,
        frames_analyzed=len(selected_frames),
        summary=summary,
        events=_state_transition_events(selected_frames, assignments),
        sampled_values=_state_samples(selected_frames, rmsd_values, spread_values, assignments),
        limitations=list(
            dict.fromkeys(
                [
                    definition("conformational_states").description,
                    *definition("conformational_states").limitations,
                    *resolved.limitations,
                ]
            )
        ),
    )


def _rmsd_timeline(
    world: InspectionWorld,
    selection: str,
    selected_frames: list[int],
) -> TimelineResult:
    universe = world.universe
    resolved = world.resolve_selection(selection, frame=selected_frames[0])
    reference_positions = np.asarray(resolved.atomgroup.positions, dtype=float).copy()
    values: list[float] = []
    for frame in selected_frames:
        universe.trajectory[frame]
        values.append(aligned_rmsd_A(reference_positions, resolved.atomgroup.positions))

    representative_frames = representative_numeric_frames(selected_frames, values, "rmsd_A")
    summary = numeric_summary(values, "rmsd_A")
    summary.update(
        {
            "reference_frame": selected_frames[0],
            "alignment": "kabsch_over_selection",
            "representative_frames": representative_frames,
        }
    )
    events = []
    high = representative_frames.get("highest")
    if high and high.get("rmsd_A", 0.0) > 0:
        events.append({"type": "rmsd_peak", **high})

    return _timeline_result(
        metric="rmsd",
        selection=selection,
        frames_analyzed=len(selected_frames),
        summary=summary,
        events=events,
        sampled_values=sampled_values(selected_frames, values, "rmsd_A"),
        limitations=list(
            dict.fromkeys(
                [
                    "RMSD is Kabsch-aligned over the selected atoms against the first selected frame.",
                    "No domain decomposition or symmetry correction is applied.",
                    *resolved.limitations,
                ]
            )
        ),
    )


def _mobility_timeline(
    world: InspectionWorld,
    metric_name: str,
    selection: str,
    selected_frames: list[int],
) -> TimelineResult:
    universe = world.universe
    resolved = world.resolve_selection(selection, frame=selected_frames[0])
    position_frames: list[np.ndarray] = []
    for frame in selected_frames:
        universe.trajectory[frame]
        position_frames.append(np.asarray(resolved.atomgroup.positions, dtype=float).copy())

    atom_rmsf = rmsf_per_atom_A(position_frames, align=True)
    residue_samples = _residue_mobility_samples(resolved.atomgroup, atom_rmsf)
    summary = numeric_summary(atom_rmsf, "rmsf_A")
    summary.update(
        {
            "reference_frame": selected_frames[0],
            "alignment": "kabsch_over_selection",
            "representative_frames": representative_displacement_frames(
                selected_frames,
                position_frames,
            ),
        }
    )
    if residue_samples:
        summary["most_mobile_object"] = residue_samples[0]["object"]
        summary["least_mobile_object"] = residue_samples[-1]["object"]

    events = []
    if residue_samples and residue_samples[0]["mean_rmsf_A"] > 0:
        events.append({"type": "mobility_peak", **residue_samples[0]})

    return _timeline_result(
        metric=metric_name,
        selection=selection,
        frames_analyzed=len(selected_frames),
        summary=summary,
        events=events,
        sampled_values=residue_samples,
        limitations=list(
            dict.fromkeys(
                [
                    "RMSF/mobility is Kabsch-aligned over the selected atoms.",
                    "Residue mobility is aggregated from selected atoms only.",
                    *resolved.limitations,
                ]
            )
        ),
    )


def _displacement_timeline(
    world: InspectionWorld,
    selection: str,
    selected_frames: list[int],
) -> TimelineResult:
    universe = world.universe
    resolved = world.resolve_selection(selection, frame=selected_frames[0])
    reference_positions = np.asarray(resolved.atomgroup.positions, dtype=float).copy()
    frame_summaries: list[dict[str, Any]] = []

    for frame in selected_frames:
        universe.trajectory[frame]
        aligned_positions = align_positions_to_reference(reference_positions, resolved.atomgroup.positions)
        atom_displacements = np.sqrt(np.sum((aligned_positions - reference_positions) ** 2, axis=1))
        residue_displacements = _residue_displacement_values(resolved.atomgroup, atom_displacements)
        most_displaced = max(
            residue_displacements,
            key=lambda item: (float(item["mean_displacement_A"]), str(item["object"])),
        )
        frame_summaries.append(
            {
                "frame": frame,
                "mean_displacement_A": round(float(np.mean(atom_displacements)), 3),
                "max_displacement_A": round(float(np.max(atom_displacements)), 3),
                "most_displaced_object": most_displaced["object"],
                "most_displaced_object_mean_A": most_displaced["mean_displacement_A"],
            }
        )

    max_values = [float(summary["max_displacement_A"]) for summary in frame_summaries]
    representative_frames = representative_numeric_frames(
        selected_frames,
        max_values,
        "max_displacement_A",
    )
    summary = numeric_summary(max_values, "displacement_A")
    summary.update(
        {
            "reference_frame": selected_frames[0],
            "alignment": "kabsch_over_selection",
            "representative_frames": representative_frames,
        }
    )
    if frame_summaries:
        peak = max(
            frame_summaries,
            key=lambda item: (float(item["max_displacement_A"]), int(item["frame"])),
        )
        summary["most_displaced_object"] = peak["most_displaced_object"]

    events = []
    high = representative_frames.get("highest")
    if high and high.get("max_displacement_A", 0.0) > 0:
        events.append({"type": "displacement_peak", **high})

    return _timeline_result(
        metric="displacement",
        selection=selection,
        frames_analyzed=len(selected_frames),
        summary=summary,
        events=events,
        sampled_values=frame_summaries,
        limitations=list(
            dict.fromkeys(
                [
                    "Displacement is Kabsch-aligned over the selected atoms against the first selected frame.",
                    "Per-residue displacement is aggregated from selected atoms only.",
                    *resolved.limitations,
                ]
            )
        ),
    )


def _selection_spread_timeline(
    world: InspectionWorld,
    selection: str,
    selected_frames: list[int],
) -> TimelineResult:
    universe = world.universe
    resolved = world.resolve_selection(selection, frame=selected_frames[0])
    spread_values: list[float] = []
    for frame in selected_frames:
        universe.trajectory[frame]
        spread_values.append(_radius_of_gyration_A(resolved.atomgroup.positions))

    representative_frames = selection_spread_representative_frames(selected_frames, spread_values)
    summary = numeric_summary(spread_values, "selection_spread_A")
    summary.update(
        {
            "spread_proxy": "selection_radius_of_gyration",
            "representative_frames": representative_frames,
        }
    )
    events = selection_spread_events(selected_frames, spread_values)

    return _timeline_result(
        metric="selection_spread",
        selection=selection,
        frames_analyzed=len(selected_frames),
        summary=summary,
        events=events,
        sampled_values=sampled_values(selected_frames, spread_values, "selection_spread_A"),
        limitations=list(
            dict.fromkeys(
                [
                    definition("selection_spread").description,
                    *definition("selection_spread").limitations,
                    *resolved.limitations,
                ]
            )
        ),
    )


def _timeline_source_id(resolved_objects: tuple[str, ...], fallback: str | None) -> str:
    if len(resolved_objects) == 1:
        return resolved_objects[0]
    if resolved_objects:
        return "selection:" + ",".join(resolved_objects[:3])
    return fallback or "selection"


def _timeline_result(**kwargs: Any) -> TimelineResult:
    summary = kwargs.get("summary", {})
    metric = str(kwargs.get("metric", "timeline"))
    frames_analyzed = int(kwargs.get("frames_analyzed", 0))
    if isinstance(summary, dict):
        kwargs.setdefault("summary_text", _timeline_summary_text(metric, summary, frames_analyzed))
    return TimelineResult(**kwargs)


def _timeline_summary_text(metric: str, summary: dict[str, Any], frames_analyzed: int) -> str:
    frame_word = "frame" if frames_analyzed == 1 else "frames"
    prefix = f"{metric} over {frames_analyzed} {frame_word}"
    if metric in {"distance", "contact", "relation"}:
        text = (
            f"{prefix}: min {_value_text(summary.get('min_A'))} A, "
            f"median {_value_text(summary.get('median_A'))} A, "
            f"max {_value_text(summary.get('max_A'))} A"
        )
        if summary.get("contact_occupancy") is not None:
            text += f"; contact occupancy {_value_text(summary.get('contact_occupancy'))}"
        if summary.get("dominant_relation_type"):
            text += f"; dominant relation {summary['dominant_relation_type']}"
        return text + "."
    if metric == "hydrogen_bonds":
        return (
            f"{prefix}: H-bond occupancy "
            f"{_value_text(summary.get('hydrogen_bond_occupancy'))}; "
            f"{summary.get('total_hydrogen_bond_observations', 0)} observations."
        )
    if metric == "interaction_persistence":
        return (
            f"{prefix}: contact occupancy {_value_text(summary.get('contact_occupancy'))}; "
            f"H-bond occupancy {_value_text(summary.get('hydrogen_bond_occupancy'))}; "
            f"dominant relation {summary.get('dominant_relation_type') or 'none'}."
        )
    if metric == "ligand_stability":
        return (
            f"{prefix}: {summary.get('stability_class') or 'unclassified'}; "
            f"max site-aligned ligand RMSD "
            f"{_value_text(summary.get('max_site_aligned_ligand_rmsd_A'))} A; "
            f"contact occupancy {_value_text(summary.get('contact_occupancy'))}."
        )
    if metric == "states":
        return (
            f"{prefix}: {summary.get('state_count', 0)} representative state(s) by "
            f"{summary.get('state_method') or 'unspecified method'}; RMSD range "
            f"{_value_text(summary.get('rmsd_range_A'))} A."
        )
    if metric == "centroid_distance":
        return (
            f"{prefix}: min {_value_text(summary.get('min_centroid_distance_A'))} A, "
            f"median {_value_text(summary.get('median_centroid_distance_A'))} A, "
            f"max {_value_text(summary.get('max_centroid_distance_A'))} A."
        )
    if metric == "rmsd":
        return (
            f"{prefix}: max RMSD {_value_text(summary.get('max_rmsd_A'))} A "
            f"against frame {summary.get('reference_frame')}."
        )
    if metric in {"rmsf", "mobility"}:
        return (
            f"{prefix}: max RMSF {_value_text(summary.get('max_rmsf_A'))} A; "
            f"most mobile object {summary.get('most_mobile_object') or 'none'}."
        )
    if metric == "displacement":
        return (
            f"{prefix}: max displacement "
            f"{_value_text(summary.get('max_displacement_A'))} A; "
            f"most displaced object {summary.get('most_displaced_object') or 'none'}."
        )
    if metric == "selection_spread":
        return (
            f"{prefix}: min {_value_text(summary.get('min_selection_spread_A'))} A, "
            f"median {_value_text(summary.get('median_selection_spread_A'))} A, "
            f"max {_value_text(summary.get('max_selection_spread_A'))} A."
        )
    return f"{prefix}: see summary, events, sampled_values, and limitations."


def _value_text(value: Any) -> str:
    if value is None:
        return "not available"
    if isinstance(value, float):
        return f"{value:g}"
    return str(value)


def _residue_mobility_samples(
    atomgroup: Any,
    atom_rmsf: list[float],
    max_samples: int = 25,
) -> list[dict[str, Any]]:
    grouped_values: dict[str, list[float]] = {}
    for atom, value in zip(atomgroup, atom_rmsf, strict=True):
        object_id = object_id_for_residue(atom.residue)
        grouped_values.setdefault(object_id, []).append(value)

    samples: list[dict[str, Any]] = []
    for object_id, values in grouped_values.items():
        numeric = np.asarray(values, dtype=float)
        samples.append(
            {
                "object": object_id,
                "mean_rmsf_A": round(float(np.mean(numeric)), 3),
                "max_atom_rmsf_A": round(float(np.max(numeric)), 3),
                "selected_atom_count": len(values),
            }
        )
    samples.sort(key=lambda item: (-float(item["mean_rmsf_A"]), str(item["object"])))
    return samples[:max_samples]


def _residue_displacement_values(atomgroup: Any, atom_displacements: np.ndarray) -> list[dict[str, Any]]:
    grouped_values: dict[str, list[float]] = {}
    for atom, value in zip(atomgroup, atom_displacements.tolist(), strict=True):
        object_id = object_id_for_residue(atom.residue)
        grouped_values.setdefault(object_id, []).append(float(value))

    values: list[dict[str, Any]] = []
    for object_id, displacements in grouped_values.items():
        numeric = np.asarray(displacements, dtype=float)
        values.append(
            {
                "object": object_id,
                "mean_displacement_A": round(float(np.mean(numeric)), 3),
                "max_atom_displacement_A": round(float(np.max(numeric)), 3),
                "selected_atom_count": len(displacements),
            }
        )
    return values


def _radius_of_gyration_A(positions: Any) -> float:
    coordinates = np.asarray(positions, dtype=float)
    center = coordinates.mean(axis=0)
    deltas = coordinates - center
    return round(float(np.sqrt(np.mean(np.sum(deltas * deltas, axis=1)))), 3)


def _rmsd_between_position_arrays(reference_positions: Any, mobile_positions: Any) -> float:
    reference = np.asarray(reference_positions, dtype=float)
    mobile = np.asarray(mobile_positions, dtype=float)
    if reference.shape != mobile.shape:
        raise MetricError("ligand stability requires consistent ligand atom correspondence.")
    deltas = mobile - reference
    return round(float(np.sqrt(np.mean(np.sum(deltas * deltas, axis=1)))), 3)


def _key_persistence_events(events: list[dict[str, Any]], max_events: int = 20) -> list[dict[str, Any]]:
    wanted = {
        "contact_forms",
        "contact_breaks",
        "relation_formed",
        "relation_lost",
        "relation_changed",
        "hydrogen_bond_forms",
        "hydrogen_bond_breaks",
        "hydrogen_bond_count_changes",
    }
    return [event for event in events if event.get("type") in wanted][:max_events]


def _interaction_persistence_samples(
    relation_samples: list[dict[str, Any]],
    hbond_samples: list[dict[str, Any]],
) -> list[dict[str, Any]]:
    hbond_by_frame = {int(sample["frame"]): sample for sample in hbond_samples}
    samples: list[dict[str, Any]] = []
    for relation_sample in relation_samples:
        frame = int(relation_sample["frame"])
        hbond_sample = hbond_by_frame.get(frame, {})
        distance = relation_sample.get("distance_A") or relation_sample.get("value_A")
        contact = relation_sample.get("contact")
        if contact is None and isinstance(distance, (int, float)):
            contact = float(distance) <= CONTACT_CUTOFF_A
        samples.append(
            {
                "frame": frame,
                "relation_type": relation_sample.get("relation_type"),
                "distance_A": distance,
                "contact": contact,
                "hydrogen_bond_count": hbond_sample.get("hydrogen_bond_count"),
                "has_hydrogen_bond": hbond_sample.get("has_hydrogen_bond"),
            }
        )
    return samples


def _ligand_stability_class(contact_occupancy: float, max_rmsd_A: float) -> str:
    if (
        contact_occupancy >= LIGAND_STABILITY_STABLE_MIN_CONTACT_OCCUPANCY
        and max_rmsd_A <= LIGAND_STABILITY_STABLE_MAX_RMSD_A
    ):
        return "stable_contact_persistent"
    if (
        contact_occupancy <= LIGAND_STABILITY_UNSTABLE_MAX_CONTACT_OCCUPANCY
        or max_rmsd_A >= LIGAND_STABILITY_UNSTABLE_MIN_RMSD_A
    ):
        return "unstable_or_repositioned"
    return "partially_stable"


def _ligand_stability_representative_frames(
    frames: list[int],
    ligand_rmsd_values: list[float],
    center_distances: list[float | None],
) -> dict[str, Any]:
    representatives = representative_numeric_frames(
        frames,
        ligand_rmsd_values,
        "site_aligned_ligand_rmsd_A",
    )
    if "lowest" in representatives:
        representatives["most_stable_pose"] = representatives.pop("lowest")
    if "highest" in representatives:
        representatives["least_stable_pose"] = representatives.pop("highest")
    distance_representatives = representative_distance_frames(frames, center_distances, [])
    if "closest" in distance_representatives:
        representatives["closest_site_centroid"] = distance_representatives["closest"]
    if "farthest" in distance_representatives:
        representatives["farthest_site_centroid"] = distance_representatives["farthest"]
    return representatives


def _ligand_stability_events(
    frames: list[int],
    ligand_rmsd_values: list[float],
    relation_events: list[dict[str, Any]],
    ligand_id: str,
    site_id: str,
) -> list[dict[str, Any]]:
    events = _key_persistence_events(relation_events)
    for frame, rmsd_value in zip(frames, ligand_rmsd_values, strict=True):
        if rmsd_value >= LIGAND_STABILITY_UNSTABLE_MIN_RMSD_A:
            events.append(
                {
                    "type": "ligand_pose_repositioned",
                    "frame": frame,
                    "site_aligned_ligand_rmsd_A": rmsd_value,
                    "ligand": ligand_id,
                    "site": site_id,
                }
            )
            break
    return events


def _ligand_stability_samples(
    frames: list[int],
    ligand_rmsd_values: list[float],
    center_distances: list[float | None],
    relation_samples: list[dict[str, Any]],
    hbond_samples: list[dict[str, Any]],
    max_samples: int = 50,
) -> list[dict[str, Any]]:
    relation_by_frame = {int(sample["frame"]): sample for sample in relation_samples}
    hbond_by_frame = {int(sample["frame"]): sample for sample in hbond_samples}
    if len(frames) <= max_samples:
        chosen = list(range(len(frames)))
    else:
        chosen = sorted(set(np.linspace(0, len(frames) - 1, max_samples, dtype=int).tolist()))
    samples: list[dict[str, Any]] = []
    for index in chosen:
        frame = frames[index]
        relation_sample = relation_by_frame.get(frame, {})
        hbond_sample = hbond_by_frame.get(frame, {})
        samples.append(
            {
                "frame": frame,
                "site_aligned_ligand_rmsd_A": ligand_rmsd_values[index],
                "ligand_site_centroid_distance_A": center_distances[index],
                "relation_type": relation_sample.get("relation_type"),
                "distance_A": relation_sample.get("distance_A"),
                "hydrogen_bond_count": hbond_sample.get("hydrogen_bond_count"),
                "has_hydrogen_bond": hbond_sample.get("has_hydrogen_bond"),
            }
        )
    return samples


def _state_assignments(rmsd_values: list[float]) -> list[str]:
    if not rmsd_values:
        return []
    rmsd_range = max(rmsd_values) - min(rmsd_values)
    if rmsd_range < STATE_MIN_RMSD_RANGE_A:
        return ["state_1"] * len(rmsd_values)
    median = float(np.median(np.asarray(rmsd_values, dtype=float)))
    assignments = ["state_1" if value <= median else "state_2" for value in rmsd_values]
    if len(set(assignments)) == 1:
        assignments[-1] = "state_2"
    return assignments


def _state_summaries(
    frames: list[int],
    rmsd_values: list[float],
    spread_values: list[float],
    assignments: list[str],
) -> list[dict[str, Any]]:
    states: list[dict[str, Any]] = []
    for state_id in sorted(set(assignments)):
        indices = [index for index, value in enumerate(assignments) if value == state_id]
        state_rmsd = [rmsd_values[index] for index in indices]
        state_spread = [spread_values[index] for index in indices]
        mean_rmsd = float(np.mean(state_rmsd))
        representative_index = min(indices, key=lambda index: abs(rmsd_values[index] - mean_rmsd))
        states.append(
            {
                "state_id": state_id,
                "frame_count": len(indices),
                "occupancy": round(len(indices) / len(frames), 3),
                "first_frame": frames[indices[0]],
                "last_frame": frames[indices[-1]],
                "representative_frame": frames[representative_index],
                "mean_rmsd_A": round(mean_rmsd, 3),
                "min_rmsd_A": round(float(np.min(state_rmsd)), 3),
                "max_rmsd_A": round(float(np.max(state_rmsd)), 3),
                "mean_selection_spread_A": round(float(np.mean(state_spread)), 3),
            }
        )
    return states


def _state_transition_events(frames: list[int], assignments: list[str]) -> list[dict[str, Any]]:
    events: list[dict[str, Any]] = []
    previous: str | None = None
    for frame, state_id in zip(frames, assignments, strict=True):
        if previous is not None and state_id != previous:
            events.append(
                {
                    "type": "state_transition",
                    "frame": frame,
                    "from_state": previous,
                    "to_state": state_id,
                }
            )
        previous = state_id
    return events


def _state_samples(
    frames: list[int],
    rmsd_values: list[float],
    spread_values: list[float],
    assignments: list[str],
    max_samples: int = 50,
) -> list[dict[str, Any]]:
    if len(frames) <= max_samples:
        chosen = list(range(len(frames)))
    else:
        chosen = sorted(set(np.linspace(0, len(frames) - 1, max_samples, dtype=int).tolist()))
    return [
        {
            "frame": frames[index],
            "state_id": assignments[index],
            "aligned_rmsd_A": rmsd_values[index],
            "selection_spread_A": spread_values[index],
        }
        for index in chosen
    ]


def _hydrogen_bond_events(
    frames: list[int],
    counts: list[int],
    source_id: str,
    target_id: str,
) -> list[dict[str, Any]]:
    events: list[dict[str, Any]] = []
    previous_present: bool | None = None
    previous_count: int | None = None
    for frame, count in zip(frames, counts, strict=True):
        present = count > 0
        if previous_present is not None and present != previous_present:
            events.append(
                {
                    "type": "hydrogen_bond_forms" if present else "hydrogen_bond_breaks",
                    "frame": frame,
                    "hydrogen_bond_count": count,
                    "source": source_id,
                    "target": target_id,
                }
            )
        elif previous_count is not None and count != previous_count and (present or previous_present):
            events.append(
                {
                    "type": "hydrogen_bond_count_changes",
                    "frame": frame,
                    "from_count": previous_count,
                    "to_count": count,
                    "source": source_id,
                    "target": target_id,
                }
            )
        previous_present = present
        previous_count = count
    return events


def _hydrogen_bond_representative_frames(
    frames: list[int],
    counts: list[int],
) -> dict[str, Any]:
    representatives = representative_numeric_frames(frames, counts, "hydrogen_bond_count")
    for frame, count in zip(frames, counts, strict=True):
        if count > 0:
            representatives["first_present"] = {"frame": frame, "hydrogen_bond_count": count}
            break
    for frame, count in zip(frames, counts, strict=True):
        if count == 0:
            representatives.setdefault("first_absent", {"frame": frame, "hydrogen_bond_count": count})
            break
    return representatives


def _sampled_hydrogen_bonds(
    universe: Any,
    frames: list[int],
    series: HydrogenBondSeries,
    max_samples: int = 50,
) -> list[dict[str, Any]]:
    if len(frames) <= max_samples:
        chosen = list(range(len(frames)))
    else:
        chosen = sorted(set(np.linspace(0, len(frames) - 1, max_samples, dtype=int).tolist()))
    grouped = series.observations_by_frame()
    samples: list[dict[str, Any]] = []
    for index in chosen:
        frame = frames[index]
        observations = sorted(
            grouped.get(frame, ()),
            key=lambda item: (-item.angle_deg, item.distance_A, item.donor_index),
        )
        sample: dict[str, Any] = {
            "frame": frame,
            "hydrogen_bond_count": len(observations),
            "has_hydrogen_bond": bool(observations),
        }
        if observations:
            best = observations[0]
            sample.update(
                {
                    "donor_atom": atom_label(universe.atoms[best.donor_index]),
                    "hydrogen_atom": atom_label(universe.atoms[best.hydrogen_index]),
                    "acceptor_atom": atom_label(universe.atoms[best.acceptor_index]),
                    "distance_A": best.distance_A,
                    "angle_deg": best.angle_deg,
                }
            )
        samples.append(sample)
    return samples


def _top_hydrogen_bond_pairs(
    universe: Any,
    series: HydrogenBondSeries,
    limit: int = 10,
) -> list[dict[str, Any]]:
    pairs = []
    for pair in series.pair_occupancies()[:limit]:
        pairs.append(
            {
                "donor_atom": atom_label(universe.atoms[int(pair["donor_index"])]),
                "acceptor_atom": atom_label(universe.atoms[int(pair["acceptor_index"])]),
                "frame_count": pair["frame_count"],
                "occupancy": pair["occupancy"],
            }
        )
    return pairs


def _relation_from_hbond_observation(
    universe: Any,
    observation: HydrogenBondObservation,
    source_atoms: Any,
    target_atoms: Any,
    source_id: str,
    target_id: str,
) -> Relation:
    hbond_definition = definition("hydrogen_bond")
    occupancy_definition = definition("hydrogen_bond_occupancy")
    source_atom_index, target_atom_index = _source_target_indices_for_hbond(
        observation,
        {int(atom.ix) for atom in source_atoms},
        {int(atom.ix) for atom in target_atoms},
    )
    donor_atom = atom_label(universe.atoms[observation.donor_index])
    hydrogen_atom = atom_label(universe.atoms[observation.hydrogen_index])
    acceptor_atom = atom_label(universe.atoms[observation.acceptor_index])
    return Relation(
        source=source_id,
        target=target_id,
        type="hydrogen_bond",
        category=hbond_definition.category,
        confidence="backend",
        backend=MDANALYSIS_HBOND_BACKEND,
        definition_id="hydrogen_bond",
        definition_source=hbond_definition.source,
        reference_keys=list(
            dict.fromkeys((*hbond_definition.reference_keys, *occupancy_definition.reference_keys))
        ),
        min_distance_A=observation.distance_A,
        cutoff_A=TEMPORAL_HBOND_DONOR_ACCEPTOR_DISTANCE_A,
        angle_deg=observation.angle_deg,
        source_atom=atom_label(universe.atoms[source_atom_index]),
        target_atom=atom_label(universe.atoms[target_atom_index]),
        evidence=[
            EvidenceItem(
                type="metric",
                metric="donor_acceptor_distance",
                value=observation.distance_A,
                unit="angstrom",
                source=f"{donor_atom}->{acceptor_atom}",
            ),
            EvidenceItem(
                type="metric",
                metric="donor_hydrogen_acceptor_angle",
                value=observation.angle_deg,
                unit="degree",
                source=f"{donor_atom}->{hydrogen_atom}->{acceptor_atom}",
            ),
            EvidenceItem(type="object", metric="donor_atom", value=donor_atom, source=donor_atom),
            EvidenceItem(
                type="object",
                metric="hydrogen_atom",
                value=hydrogen_atom,
                source=hydrogen_atom,
            ),
            EvidenceItem(
                type="object",
                metric="acceptor_atom",
                value=acceptor_atom,
                source=acceptor_atom,
            ),
        ],
        limitations=list(occupancy_definition.limitations),
    )


def _source_target_indices_for_hbond(
    observation: HydrogenBondObservation,
    source_indices: set[int],
    target_indices: set[int],
) -> tuple[int, int]:
    if observation.donor_index in source_indices and observation.acceptor_index in target_indices:
        return observation.donor_index, observation.acceptor_index
    if observation.acceptor_index in source_indices and observation.donor_index in target_indices:
        return observation.acceptor_index, observation.donor_index
    if observation.donor_index in source_indices:
        return observation.donor_index, observation.acceptor_index
    if observation.acceptor_index in source_indices:
        return observation.acceptor_index, observation.donor_index
    return observation.donor_index, observation.acceptor_index


def _relation_priority(relation_type: str) -> int:
    return relation_priority_for_type(relation_type)


def _sampled_relations(
    frames: list[int],
    relations: list[Relation | None],
    max_samples: int = 50,
) -> list[dict[str, Any]]:
    if len(frames) <= max_samples:
        chosen = list(range(len(frames)))
    else:
        step = max(1, (len(frames) - 1) // (max_samples - 1))
        chosen = list(range(0, len(frames), step))[: max_samples - 1]
        if chosen[-1] != len(frames) - 1:
            chosen.append(len(frames) - 1)

    samples: list[dict[str, Any]] = []
    for index in chosen:
        relation = relations[index]
        if relation is None:
            samples.append({"frame": frames[index], "relation_type": None})
            continue
        samples.append(
            {
                "frame": frames[index],
                "relation_type": relation.type,
                "distance_A": relation.min_distance_A,
                "source_atom": relation.source_atom,
                "target_atom": relation.target_atom,
            }
        )
    return samples
