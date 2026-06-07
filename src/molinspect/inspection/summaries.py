"""Summary text and small aggregate helpers for inspection outputs."""

from __future__ import annotations

from typing import Any

from ..definitions import CONTACT_CUTOFF_A
from ..schemas import Relation, StructuralProfile


def location_summary(resolved_objects: tuple[str, ...], location: dict[str, Any]) -> str:
    """Return a compact human-readable location summary from structured fields."""

    label = ", ".join(resolved_objects[:3]) if resolved_objects else "selection"
    if len(resolved_objects) > 3:
        label += f", and {len(resolved_objects) - 3} more"
    chain = location.get("chain")
    nearest = location.get("nearest_protein_residues") or []
    ligands = location.get("near_ligands") or []
    profile = location.get("structural_profile")

    chain_text = format_chain_text(chain)
    parts = [f"{label} is in {chain_text}." if chain_text else f"{label} was resolved."]
    if isinstance(profile, StructuralProfile) and profile.secondary_structure:
        parts.append(f"Secondary structure is {profile.secondary_structure}.")
    if isinstance(profile, StructuralProfile) and profile.exposure:
        source = profile.exposure_source
        source_text = f" ({source})" if source else ""
        parts.append(f"Exposure is {profile.exposure}{source_text}.")
    if isinstance(profile, StructuralProfile) and profile.local_packing:
        parts.append(f"Local packing is {profile.local_packing}.")
    if isinstance(profile, StructuralProfile) and profile.interface_chains:
        chains = ", ".join(profile.interface_chains)
        parts.append(f"Interface contacts involve chain(s) {chains}.")
    if isinstance(profile, StructuralProfile) and profile.biological_interface_ids:
        parts.append("PISA biological interface membership is present.")
    if nearest:
        parts.append(
            f"Nearest protein residue is {nearest[0].object_id} at "
            f"{nearest[0].min_distance_A} A."
        )
    if ligands:
        parts.append(
            f"Nearest ligand/ion is {ligands[0].object_id} at "
            f"{ligands[0].min_distance_A} A."
        )
    return " ".join(parts)


def relation_counts(relations: list[Relation]) -> dict[str, int]:
    """Return sorted relation-type counts."""

    counts: dict[str, int] = {}
    for relation in relations:
        counts[relation.type] = counts.get(relation.type, 0) + 1
    return dict(sorted(counts.items()))


def context_summary(
    relations: list[Relation],
    radius: float,
    contact_count: int,
    relation_type_counts: dict[str, int],
) -> str:
    """Return a compact summary of nearby context relations."""

    if not relations:
        return f"No nearby non-water objects within {radius:g} A."
    relation_text = ", ".join(f"{count} {name}" for name, count in relation_type_counts.items())
    closest = min(
        relations,
        key=lambda relation: relation.min_distance_A
        if relation.min_distance_A is not None
        else float("inf"),
    )
    top_relation = relations[0]
    top_distance = (
        f" at {top_relation.min_distance_A} A" if top_relation.min_distance_A is not None else ""
    )
    contact_phrase = "1 is a contact" if contact_count == 1 else f"{contact_count} are contacts"
    return (
        f"{len(relations)} nearby objects within {radius:g} A; "
        f"{contact_phrase} at <= {CONTACT_CUTOFF_A:g} A. "
        f"Relation types: {relation_text}. "
        f"Top-ranked relation is {top_relation.type} to {top_relation.target}{top_distance}. "
        f"Closest target is {closest.target} at {closest.min_distance_A} A."
    )


def format_chain_text(chain: Any) -> str:
    """Format one or more chain IDs for location summaries."""

    if isinstance(chain, list):
        if not chain:
            return ""
        if len(chain) == 1:
            return f"chain {chain[0]}"
        return "chains " + ", ".join(str(value) for value in chain)
    if chain:
        return f"chain {chain}"
    return ""
