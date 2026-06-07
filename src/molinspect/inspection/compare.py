"""Compare-state inspection primitive."""

from __future__ import annotations

from typing import Any

import numpy as np

from ..definitions import CONTACT_CUTOFF_A, DEFAULT_CONTEXT_RADIUS_A, DISTANCE_CHANGE_MIN_DELTA_A
from ..errors import MetricError
from ..metrics import aligned_rmsd_A
from ..schemas import CompareResult, EvidenceItem, ObjectRef
from ..world import InspectionWorld
from .events import relation_change_type
from .static import context


def compare(
    world: InspectionWorld,
    selection: str,
    frame_a: int | str,
    frame_b: int | str,
    radius: float = DEFAULT_CONTEXT_RADIUS_A,
) -> CompareResult:
    """Compare compact local context around a selection between two frames."""

    context_a = context(world, selection, frame=frame_a, radius=radius)
    context_b = context(world, selection, frame=frame_b, radius=radius)

    rels_a = {relation.target: relation for relation in context_a.relations}
    rels_b = {relation.target: relation for relation in context_b.relations}
    objects_a = {obj.id: obj for obj in context_a.objects}
    objects_b = {obj.id: obj for obj in context_b.objects}

    targets_a = set(rels_a)
    targets_b = set(rels_b)
    gained = sorted(targets_b - targets_a)
    lost = sorted(targets_a - targets_b)
    common = sorted(targets_a & targets_b)

    changes: list[dict[str, Any]] = []
    for object_id in gained:
        changes.append(
            {
                "type": _nearby_gain_type(selection, objects_b.get(object_id)),
                "object": object_id,
                "frame": context_b.frame,
            }
        )
    for object_id in lost:
        changes.append(
            {
                "type": _nearby_loss_type(selection, objects_a.get(object_id)),
                "object": object_id,
                "frame": context_b.frame,
            }
        )

    for object_id in common:
        rel_a = rels_a[object_id]
        rel_b = rels_b[object_id]
        if rel_a.min_distance_A is None or rel_b.min_distance_A is None:
            continue
        delta = round(rel_b.min_distance_A - rel_a.min_distance_A, 3)
        if abs(delta) >= DISTANCE_CHANGE_MIN_DELTA_A:
            changes.append(
                {
                    "type": "distance_change",
                    "object": object_id,
                    "from_A": rel_a.min_distance_A,
                    "to_A": rel_b.min_distance_A,
                    "delta_A": delta,
                }
            )
        if rel_a.type != rel_b.type:
            changes.append(
                {
                    "type": relation_change_type(rel_a.type, rel_b.type),
                    "object": object_id,
                    "from_relation": rel_a.type,
                    "to_relation": rel_b.type,
                }
            )

    rmsd_A = _selection_rmsd_between_frames(world, selection, context_a.frame, context_b.frame)
    evidence = [
        EvidenceItem(type="parameter", metric="radius", value=radius, unit="angstrom"),
        EvidenceItem(
            type="parameter",
            metric="contact_cutoff",
            value=CONTACT_CUTOFF_A,
            unit="angstrom",
            source="nonbonded_contact",
        ),
    ]
    if rmsd_A is not None:
        changes.append(
            {
                "type": "selection_rmsd",
                "from_frame": context_a.frame,
                "to_frame": context_b.frame,
                "rmsd_A": rmsd_A,
                "alignment": "kabsch_over_selection",
            }
        )
        evidence.append(
            EvidenceItem(
                type="metric",
                metric="selection_rmsd",
                value=rmsd_A,
                unit="angstrom",
            )
        )

    changes.sort(key=_change_priority)
    return CompareResult(
        selection=selection,
        frame_a=context_a.frame,
        frame_b=context_b.frame,
        main_changes=changes[:25],
        summary=(
            f"{len(gained)} nearby objects gained, {len(lost)} lost, "
            f"{len(common)} shared within {radius:g} A"
            + (f"; selection RMSD {rmsd_A} A." if rmsd_A is not None else ".")
        ),
        evidence=evidence,
        limitations=sorted(set(context_a.limitations + context_b.limitations)),
    )


def _nearby_gain_type(selection: str, obj: ObjectRef | None) -> str:
    if _selection_or_object_is_interface(selection, obj):
        return "interface_contact_gain"
    return "nearby_object_gain"


def _nearby_loss_type(selection: str, obj: ObjectRef | None) -> str:
    if _selection_or_object_is_interface(selection, obj):
        return "interface_contact_loss"
    return "nearby_object_loss"


def _selection_or_object_is_interface(selection: str, obj: ObjectRef | None) -> bool:
    if selection.strip().startswith(("interchain_contact_interface:", "biological_interface:")):
        return True
    if obj is None:
        return False
    return bool(
        obj.annotations.get("interchain_contact_interface_ids")
        or obj.annotations.get("biological_interface_ids")
        or obj.id.startswith("interchain_contact_interface:")
        or obj.id.startswith("biological_interface:")
    )


def _selection_rmsd_between_frames(
    world: InspectionWorld,
    selection: str,
    frame_a: int,
    frame_b: int,
) -> float | None:
    try:
        resolved_a = world.resolve_selection(selection, frame=frame_a)
        positions_a = np.asarray(resolved_a.atomgroup.positions, dtype=float).copy()
        resolved_b = world.resolve_selection(selection, frame=frame_b)
        return aligned_rmsd_A(positions_a, resolved_b.atomgroup.positions)
    except MetricError:
        return None


def _change_priority(change: dict[str, Any]) -> tuple[int, float, str]:
    type_priority = {
        "contact_formed": 0,
        "contact_lost": 0,
        "interface_contact_gain": 0,
        "interface_contact_loss": 0,
        "ligand_moved_toward": 0,
        "ligand_moved_away": 0,
        "nearby_object_gain": 1,
        "nearby_object_loss": 1,
        "distance_change": 2,
        "selection_rmsd": 3,
    }.get(str(change.get("type")), 9)
    delta = abs(float(change.get("delta_A", 0.0)))
    return (type_priority, -delta, str(change.get("object", "")))
