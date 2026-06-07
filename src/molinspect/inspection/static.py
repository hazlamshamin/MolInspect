"""Static spatial inspection primitives."""

from __future__ import annotations

from collections.abc import Sequence
from typing import Any

from ..definitions import CONTACT_CUTOFF_A
from ..errors import MetricError
from ..metrics import center_of_geometry_A, distance_between_centers_A
from ..objects import chain_identifier_for_atom
from ..schemas import (
    ContextBudget,
    ContextFocus,
    ContextResult,
    ContextScale,
    EvidenceItem,
    LocationInfo,
    LocateResult,
    Relation,
    StructuralProfile,
)
from ..world import InspectionWorld, normalize_context_focus
from .scales import resolve_context_settings
from .summaries import context_summary, location_summary, relation_counts
from .validation import budget_limit


def locate(
    world: InspectionWorld,
    selection: str,
    frame: int | str = 0,
    include_metrics: bool = True,
) -> LocateResult:
    """Answer where a selected object or region is, using computed evidence only."""

    universe = world.universe
    resolved = world.resolve_selection(selection, frame=frame)
    selected = resolved.atomgroup
    selected_refs = world.object_refs_for_atomgroup(selected)
    chains = _chains_for_atomgroup(selected)
    segments = _segments_for_atomgroup(selected)
    residues = [ref.id for ref in selected_refs if ref.type != "atom"]

    location: dict[str, Any] = {
        "center_of_geometry_A": center_of_geometry_A(selected),
        "selected_atom_count": len(selected),
        "selected_object_count": len(selected_refs),
        "chain": chains[0] if len(chains) == 1 else chains,
        "segment": segments[0] if len(segments) == 1 else segments,
        "objects": residues or resolved.resolved_objects,
    }
    structural_profile = world.selected_structural_profile(selected)
    if structural_profile:
        location["structural_profile"] = structural_profile
        _promote_single_profile_fields(location, structural_profile)

    evidence: list[EvidenceItem] = [
        EvidenceItem(
            type="metric",
            metric="center_of_geometry",
            value=location["center_of_geometry_A"],
            unit="angstrom",
            frame=resolved.frame,
        )
    ]

    if include_metrics:
        nearest_residues = world.nearest_entries(selected, object_types={"residue"}, limit=5)
        nearest_ligands = world.nearest_entries(
            selected,
            object_types={"ligand", "ion"},
            limit=5,
        )
        location["nearest_protein_residues"] = nearest_residues
        location["near_ligands"] = nearest_ligands
        chain_distance = _distance_to_chain_centroid_A(universe, selected, chains)
        location["distance_to_chain_centroid_A"] = chain_distance

        if nearest_residues:
            evidence.append(
                EvidenceItem(
                    type="metric",
                    metric="nearest_protein_residue_distance",
                    value=nearest_residues[0].min_distance_A,
                    unit="angstrom",
                    frame=resolved.frame,
                    source=nearest_residues[0].object_id,
                )
            )
        if nearest_ligands:
            evidence.append(
                EvidenceItem(
                    type="metric",
                    metric="nearest_ligand_distance",
                    value=nearest_ligands[0].min_distance_A,
                    unit="angstrom",
                    frame=resolved.frame,
                    source=nearest_ligands[0].object_id,
                )
            )
        if chain_distance is not None:
            evidence.append(
                EvidenceItem(
                    type="metric",
                    metric="distance_to_chain_centroid",
                    value=chain_distance,
                    unit="angstrom",
                    frame=resolved.frame,
                )
            )

    location["plain_language"] = location_summary(resolved.resolved_objects, location)

    limitations = list(resolved.limitations)
    limitations.extend(world.annotations.limitations)
    if structural_profile is None or not structural_profile.ligand_contacts:
        limitations.append(
            "No ligand-contact-shell membership was found for this selection; "
            "pocket objects require a dedicated pocket backend."
        )

    return LocateResult(
        selection=selection,
        resolved_objects=list(resolved.resolved_objects),
        frame=resolved.frame,
        location=LocationInfo(**location),
        evidence=evidence,
        limitations=limitations,
    )


def context(
    world: InspectionWorld,
    selection: str,
    frame: int | str = 0,
    radius: float | None = None,
    budget: ContextBudget | str | None = None,
    focus: ContextFocus | Sequence[ContextFocus] | None = None,
    scale: ContextScale | str | None = None,
) -> ContextResult:
    """Retrieve compact radius-based context around a selection."""

    settings = resolve_context_settings(scale=scale, radius=radius, budget=budget, focus=focus)
    if settings.radius_A <= 0:
        raise MetricError("radius must be > 0")
    limit = budget_limit(settings.budget)
    context_focus = normalize_context_focus(settings.focus)
    resolved = world.resolve_selection(selection, frame=frame)
    selected_refs = world.object_refs_for_atomgroup(resolved.atomgroup)
    selected_ids = {ref.id for ref in selected_refs}

    nearby, truncated = world.nearby_entries(
        resolved.atomgroup,
        radius_A=settings.radius_A,
        limit=limit,
        focus=context_focus,
    )

    objects_by_id = {ref.id: ref for ref in selected_refs}
    relations: list[Relation] = []
    source = selected_refs[0].id if selected_refs else "selection"
    for entry in nearby:
        ref = entry.object
        if ref.id in selected_ids:
            continue
        objects_by_id.setdefault(ref.id, ref)
        relations.append(entry.relation.model_copy(update={"source": source}))

    contact_count = sum(
        1
        for relation in relations
        if relation.min_distance_A is not None and relation.min_distance_A <= CONTACT_CUTOFF_A
    )
    relation_type_counts = relation_counts(relations)
    limitations = list(resolved.limitations)
    limitations.extend(world.annotations.limitations)
    limitations.extend(world.interaction_backend_limitations)
    if truncated:
        limitations.append(
            f"Context truncated to budget={settings.budget!r} ({limit} nearby objects)."
        )

    evidence = [
        EvidenceItem(
            type="parameter",
            metric="radius",
            value=settings.radius_A,
            unit="angstrom",
            frame=resolved.frame,
        ),
        EvidenceItem(
            type="parameter",
            metric="budget",
            value=settings.budget,
            frame=resolved.frame,
        ),
        EvidenceItem(
            type="parameter",
            metric="contact_cutoff",
            value=CONTACT_CUTOFF_A,
            unit="angstrom",
            frame=resolved.frame,
            source="nonbonded_contact",
        ),
    ]
    if settings.scale is not None:
        evidence.insert(
            0,
            EvidenceItem(
                type="parameter",
                metric="context_scale",
                value=settings.scale,
                frame=resolved.frame,
            ),
        )

    return ContextResult(
        selection=selection,
        frame=resolved.frame,
        scale=settings.scale,
        radius_A=settings.radius_A,
        budget=settings.budget,
        focus=list(context_focus),
        objects=list(objects_by_id.values()),
        relations=relations,
        summary=context_summary(relations, settings.radius_A, contact_count, relation_type_counts),
        evidence=evidence,
        limitations=limitations,
    )


def _distance_to_chain_centroid_A(universe: Any, selected: Any, chains: list[str]) -> float | None:
    if len(chains) != 1:
        return None
    chain = chains[0]
    try:
        chain_atoms = universe.select_atoms(f"chainID {chain}")
    except Exception:
        chain_atoms = universe.atoms[:0]
    if len(chain_atoms) == 0:
        try:
            chain_atoms = universe.select_atoms(f"segid {chain}")
        except Exception:
            return None
    return distance_between_centers_A(selected, chain_atoms)


def _promote_single_profile_fields(location: dict[str, Any], profile: StructuralProfile) -> None:
    """Copy common single-residue annotation fields to stable top-level keys."""

    location["secondary_structure"] = profile.secondary_structure
    location["exposure_status"] = profile.exposure
    location["exposure_source"] = profile.exposure_source
    location["surface_status"] = profile.surface_status
    location["sasa_A2"] = profile.sasa_A2
    location["relative_sasa"] = profile.relative_sasa
    location["local_packing"] = profile.local_packing
    location["local_packing_source"] = profile.local_packing_source
    location["local_contact_count"] = profile.local_contact_count
    location["local_contact_radius_A"] = profile.local_contact_radius_A
    location["interface_chains"] = profile.interface_chains
    location["nearest_interchain_distance_A"] = profile.nearest_interchain_distance_A
    location["ligand_contacts"] = profile.ligand_contacts
    location["landmark_memberships"] = profile.landmark_memberships
    location["secondary_structure_element"] = profile.secondary_structure_element
    location["ligand_contact_shell_ids"] = profile.ligand_contact_shell_ids
    location["pocket_ids"] = profile.pocket_ids
    location["interchain_contact_interface_ids"] = profile.interchain_contact_interface_ids
    location["biological_interface_ids"] = profile.biological_interface_ids


def _chains_for_atomgroup(atomgroup: Any) -> list[str]:
    return sorted({chain_identifier_for_atom(atom) for atom in atomgroup if chain_identifier_for_atom(atom)})


def _segments_for_atomgroup(atomgroup: Any) -> list[str]:
    segments = {str(getattr(atom, "segid", "")).strip() for atom in atomgroup}
    return sorted(segment for segment in segments if segment)
