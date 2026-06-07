"""Timeline event and representative-frame helpers."""

from __future__ import annotations

from collections.abc import Sequence
from typing import Any

import numpy as np

from ..definitions import (
    DISTANCE_CHANGE_MIN_DELTA_A,
    LIGAND_MOTION_MIN_DELTA_A,
    SELECTION_SPREAD_MIN_DELTA_A,
)
from ..metrics import aligned_rmsd_A
from ..schemas import Relation


def numeric_summary(values: Sequence[float], value_key: str) -> dict[str, Any]:
    """Return min/median/max/mean values for one numeric timeline metric."""

    if not values:
        return {
            f"min_{value_key}": None,
            f"median_{value_key}": None,
            f"max_{value_key}": None,
            f"mean_{value_key}": None,
        }
    numeric = np.asarray(values, dtype=float)
    return {
        f"min_{value_key}": round(float(np.min(numeric)), 3),
        f"median_{value_key}": round(float(np.median(numeric)), 3),
        f"max_{value_key}": round(float(np.max(numeric)), 3),
        f"mean_{value_key}": round(float(np.mean(numeric)), 3),
    }


def representative_numeric_frames(
    frames: list[int],
    values: Sequence[float],
    value_key: str,
) -> dict[str, dict[str, Any]]:
    """Return lowest/highest representative frames for numeric values."""

    if not values:
        return {}
    numeric = np.asarray(values, dtype=float)
    low_index = int(np.argmin(numeric))
    high_index = int(np.argmax(numeric))
    return {
        "lowest": {"frame": frames[low_index], value_key: round(float(numeric[low_index]), 3)},
        "highest": {"frame": frames[high_index], value_key: round(float(numeric[high_index]), 3)},
    }


def representative_displacement_frames(
    frames: list[int],
    position_frames: list[np.ndarray],
) -> dict[str, dict[str, Any]]:
    """Return representative aligned-RMSD frames for displacement-style timelines."""

    if not position_frames:
        return {}
    reference = position_frames[0]
    values = [aligned_rmsd_A(reference, positions) for positions in position_frames]
    return representative_numeric_frames(frames, values, "selection_rmsd_A")


def representative_distance_frames(
    frames: list[int],
    distances: list[float | None],
    events: list[dict[str, Any]],
) -> dict[str, Any]:
    """Return closest/farthest and contact-transition representative frames."""

    valid_indices = [index for index, distance in enumerate(distances) if distance is not None]
    if not valid_indices:
        return {}
    closest_index = min(valid_indices, key=lambda index: distances[index] or float("inf"))
    farthest_index = max(valid_indices, key=lambda index: distances[index] or float("-inf"))
    contact_forms = [event for event in events if event.get("type") == "contact_forms"]
    contact_breaks = [event for event in events if event.get("type") == "contact_breaks"]
    representatives: dict[str, Any] = {
        "closest": {"frame": frames[closest_index], "distance_A": distances[closest_index]},
        "farthest": {"frame": frames[farthest_index], "distance_A": distances[farthest_index]},
    }
    if contact_forms:
        representatives["contact_formed"] = contact_forms[0]
    if contact_breaks:
        representatives["contact_broken"] = contact_breaks[0]
    return representatives


def representative_relation_frames(
    frames: list[int],
    relations: list[Relation | None],
    distances: list[float | None],
) -> dict[str, Any]:
    """Return representative frames grouped by relation type."""

    valid_indices = [index for index, distance in enumerate(distances) if distance is not None]
    if not valid_indices:
        return {}
    closest_index = min(valid_indices, key=lambda index: distances[index] or float("inf"))
    farthest_index = max(valid_indices, key=lambda index: distances[index] or float("-inf"))
    by_relation_type: dict[str, int] = {}
    for index in valid_indices:
        relation = relations[index]
        if relation is not None:
            by_relation_type.setdefault(relation.type, frames[index])
    closest_relation = relations[closest_index]
    farthest_relation = relations[farthest_index]
    return {
        "closest": {
            "frame": frames[closest_index],
            "distance_A": distances[closest_index],
            "relation_type": closest_relation.type if closest_relation else None,
        },
        "farthest": {
            "frame": frames[farthest_index],
            "distance_A": distances[farthest_index],
            "relation_type": farthest_relation.type if farthest_relation else None,
        },
        "by_relation_type": by_relation_type,
    }


def movement_events(
    frames: list[int],
    distances: list[float | None],
    source_id: str,
    target_id: str,
) -> list[dict[str, Any]]:
    """Return coarse ligand moved-toward/away events from distance changes."""

    if not _involves_ligand(source_id, target_id):
        return []

    events: list[dict[str, Any]] = []
    previous_distance: float | None = None
    for frame, distance in zip(frames, distances, strict=True):
        if previous_distance is not None and distance is not None:
            delta = round(distance - previous_distance, 3)
            if abs(delta) >= LIGAND_MOTION_MIN_DELTA_A:
                events.append(
                    {
                        "type": "ligand_moved_away" if delta > 0 else "ligand_moved_toward",
                        "frame": frame,
                        "from_distance_A": previous_distance,
                        "to_distance_A": distance,
                        "delta_A": delta,
                    }
                )
        if distance is not None:
            previous_distance = distance
    return events


def centroid_distance_events(
    frames: list[int],
    distances: list[float | None],
) -> list[dict[str, Any]]:
    """Return coarse center-distance increase/decrease events."""

    events: list[dict[str, Any]] = []
    previous_distance: float | None = None
    for frame, distance in zip(frames, distances, strict=True):
        if previous_distance is not None and distance is not None:
            delta = round(distance - previous_distance, 3)
            if abs(delta) >= DISTANCE_CHANGE_MIN_DELTA_A:
                events.append(
                    {
                        "type": (
                            "centroid_distance_increased"
                            if delta > 0
                            else "centroid_distance_decreased"
                        ),
                        "frame": frame,
                        "from_centroid_distance_A": previous_distance,
                        "to_centroid_distance_A": distance,
                        "delta_A": delta,
                    }
                )
        if distance is not None:
            previous_distance = distance
    return events


def relation_change_events(
    frames: list[int],
    relations: list[Relation | None],
) -> list[dict[str, Any]]:
    """Return relation-type transition events over frames."""

    events: list[dict[str, Any]] = []
    previous: Relation | None = None
    for frame, relation in zip(frames, relations, strict=True):
        if previous is not None and relation is not None and previous.type != relation.type:
            events.append(
                {
                    "type": relation_change_type(previous.type, relation.type),
                    "frame": frame,
                    "from_relation": previous.type,
                    "to_relation": relation.type,
                    "distance_A": relation.min_distance_A,
                }
            )
        previous = relation
    return events


def relation_change_type(from_relation: str, to_relation: str) -> str:
    """Name a relation-type transition."""

    if from_relation != "nonbonded_contact" and to_relation == "nonbonded_contact":
        return "contact_formed"
    if from_relation == "nonbonded_contact" and to_relation != "nonbonded_contact":
        return "contact_lost"
    if from_relation == "near" and to_relation != "near":
        return "relation_formed"
    if from_relation != "near" and to_relation == "near":
        return "relation_lost"
    return "relation_changed"


def selection_spread_representative_frames(
    frames: list[int],
    spread_values: list[float],
) -> dict[str, dict[str, Any]]:
    """Return least/most-spread representative frames."""

    representatives = representative_numeric_frames(frames, spread_values, "selection_spread_A")
    if "lowest" in representatives:
        representatives["least_spread"] = representatives.pop("lowest")
    if "highest" in representatives:
        representatives["most_spread"] = representatives.pop("highest")
    return representatives


def selection_spread_events(
    frames: list[int],
    spread_values: list[float],
) -> list[dict[str, Any]]:
    """Return coarse selection spread increase/decrease events."""

    events: list[dict[str, Any]] = []
    previous: float | None = None
    for frame, value in zip(frames, spread_values, strict=True):
        if previous is not None:
            delta = round(value - previous, 3)
            if abs(delta) >= SELECTION_SPREAD_MIN_DELTA_A:
                events.append(
                    {
                        "type": "selection_spread_increased"
                        if delta > 0
                        else "selection_spread_decreased",
                        "frame": frame,
                        "from_selection_spread_A": previous,
                        "to_selection_spread_A": value,
                        "delta_A": delta,
                    }
                )
        previous = value
    return events


def _involves_ligand(source_id: str, target_id: str) -> bool:
    ids = (source_id, target_id)
    return any(
        object_id.startswith(("ligand:", "ion:", "ligand_contact_shell:", "pocket:"))
        or "ligand:" in object_id
        for object_id in ids
    )
