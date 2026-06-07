"""Session-owned molecular inspection world.

The public API stays intentionally small, but this module gives the internals a
single place to run structural inspections over objects, annotations, relations,
and frames.
"""

from __future__ import annotations

from collections.abc import Sequence
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, cast

from .annotations import AnnotationStore
from .errors import MetricError
from .backends.interactions import (
    INTERACTION_BACKEND_FOCUS,
    InteractionStore,
    backend_relation_or_fallback,
    build_interaction_store,
)
from .landmarks import LandmarkIndex, build_landmark_index
from .metrics import AtomPairDistance, closest_heavy_atom_pair, heavy_atomgroup, min_distance_A
from .notation import CONTEXT_FOCUS_NOTATION_HELP
from .objects import (
    chain_identifier_for_atom,
    object_refs_for_atomgroup as base_object_refs_for_atomgroup,
    object_ref_for_residue,
    residue_object_type,
)
from .relations import atom_element, chain_for_residue, relation_for_atomgroups
from .relations import relation_priority_for_type, water_bridge_relation
from .schemas import ContextFocus, NearbyObject, ObjectRef, Relation, StructuralProfile
from .selections import ResolvedSelection, resolve_selection
from .definitions import (
    HBOND_DONOR_ACCEPTOR_ELEMENTS,
    INTERFACE_DISTANCE_A,
    LIGAND_CONTACT_DISTANCE_A,
    MAX_RESOLVED_OBJECT_REFS,
    WATER_BRIDGE_MAX_DISTANCE_A,
    WATER_BRIDGE_MIN_DISTANCE_A,
)

CONTEXT_FOCUS_VALUES: tuple[ContextFocus, ...] = (
    "general",
    "contacts",
    "ligand_contact_shell",
    "metal_coordination",
    "interchain_interfaces",
    "hydrogen_bonds",
    "salt_bridges",
    "hydrophobic_contacts",
    "pi_stacking",
    "water_bridges",
    "steric_clashes",
)
DEFAULT_CONTEXT_FOCUS: tuple[ContextFocus, ...] = ("general",)
CONTEXT_RELATION_PRIORITY = {
    "metal_coordination": 0,
    "salt_bridge": 1,
    "hydrogen_bond": 2,
    "polar_contact_candidate": 3,
    "pi_stacking": 4,
    "hydrophobic_contact": 5,
    "water_bridge_candidate": 6,
    "steric_clash": 7,
    "nonbonded_contact": 8,
    "topology_bond": 9,
    "inferred_covalent_bond": 9,
    "near": 10,
}


@dataclass(slots=True)
class NearbyEntry:
    """One nearby residue-like object with structural relation evidence."""

    object: ObjectRef
    relation: Relation
    distance_A: float
    is_interchain: bool
    is_ligand_or_ion: bool
    annotation: dict[str, Any] = field(default_factory=dict)


@dataclass(slots=True)
class InspectionWorld:
    """Internal world model for one loaded structure or trajectory."""

    universe: Any
    source_files: tuple[Path, ...]
    _annotations: AnnotationStore | None = None
    _landmarks: LandmarkIndex | None = None
    _landmarks_with_pockets: LandmarkIndex | None = None
    _landmarks_with_interfaces: LandmarkIndex | None = None
    _landmarks_with_pockets_and_interfaces: LandmarkIndex | None = None
    _interaction_store: InteractionStore | None = None
    _interaction_backend_attempted: bool = False

    @property
    def annotations(self) -> AnnotationStore:
        """Return lazily computed structural annotations."""

        if self._annotations is None:
            self._annotations = AnnotationStore.build(self.universe, self.source_files)
        return self._annotations

    @property
    def landmarks(self) -> LandmarkIndex:
        """Return lazily computed biological landmark objects."""

        if self._landmarks is None:
            self._landmarks = build_landmark_index(self.universe, self.annotations)
        return self._landmarks

    def _landmark_index(
        self,
        include_backend_pockets: bool = False,
        include_backend_interfaces: bool = False,
    ) -> LandmarkIndex:
        """Return baseline landmarks, optionally enriched with backend landmarks."""

        if not include_backend_pockets and not include_backend_interfaces:
            return self.landmarks
        if include_backend_pockets and include_backend_interfaces:
            if self._landmarks_with_pockets_and_interfaces is None:
                self._landmarks_with_pockets_and_interfaces = build_landmark_index(
                    self.universe,
                    self.annotations,
                    include_backend_pockets=True,
                    include_backend_interfaces=True,
                    source_files=self.source_files,
                )
            return self._landmarks_with_pockets_and_interfaces
        if include_backend_pockets:
            if self._landmarks_with_pockets is None:
                self._landmarks_with_pockets = build_landmark_index(
                    self.universe,
                    self.annotations,
                    include_backend_pockets=True,
                    source_files=self.source_files,
                )
            return self._landmarks_with_pockets
        if self._landmarks_with_interfaces is None:
            self._landmarks_with_interfaces = build_landmark_index(
                self.universe,
                self.annotations,
                include_backend_interfaces=True,
                source_files=self.source_files,
            )
        return self._landmarks_with_interfaces

    def _active_landmarks(self) -> LandmarkIndex:
        """Return the richest requested landmark index."""

        return (
            self._landmarks_with_pockets_and_interfaces
            or self._landmarks_with_interfaces
            or self._landmarks_with_pockets
            or self.landmarks
        )

    @property
    def interaction_store(self) -> InteractionStore:
        """Return lazily computed backend interaction contacts."""

        self._interaction_backend_attempted = True
        if self._interaction_store is None:
            self._interaction_store = build_interaction_store(self.source_files)
        return self._interaction_store

    @property
    def interaction_backend_limitations(self) -> tuple[str, ...]:
        """Return interaction-backend limitations only after an attempted lookup."""

        if not self._interaction_backend_attempted or self._interaction_store is None:
            return ()
        return self._interaction_store.limitations

    def resolve_selection(self, selection: str, frame: int | str = 0) -> ResolvedSelection:
        """Resolve selections, including returned landmark object IDs."""

        landmarks = self._landmark_index(
            include_backend_pockets=selection.startswith("pocket:"),
            include_backend_interfaces=selection.startswith("biological_interface:"),
        )
        translated = landmarks.selection_by_id.get(selection, selection)
        return resolve_selection(self.universe, translated, frame=frame)

    def landmark_objects(
        self,
        include_backend_pockets: bool = False,
        include_backend_interfaces: bool = False,
    ) -> list[ObjectRef]:
        """Return computed biological landmark objects for `objects(type=...)`."""

        return list(
            self._landmark_index(
                include_backend_pockets=include_backend_pockets,
                include_backend_interfaces=include_backend_interfaces,
            ).objects
        )

    def landmark_limitations(
        self,
        include_backend_pockets: bool = False,
        include_backend_interfaces: bool = False,
    ) -> tuple[str, ...]:
        """Return limitations from requested backend landmark detectors."""

        return self._landmark_index(
            include_backend_pockets=include_backend_pockets,
            include_backend_interfaces=include_backend_interfaces,
        ).limitations

    def object_ref_for_residue(self, residue: Any, include_annotations: bool = True) -> ObjectRef:
        """Return a residue object reference with optional structural annotations."""

        if include_annotations:
            ref = self.annotations.object_ref_with_annotations(residue)
            updates = self._active_landmarks().profile_update_for_residue(residue.ix)
            annotations = {**ref.annotations, **_non_empty_landmark_updates(updates)}
            return ref.model_copy(update={"annotations": annotations})
        return object_ref_for_residue(residue)

    def object_refs_for_atomgroup(self, atomgroup: Any) -> list[ObjectRef]:
        """Return selected object refs with residue annotations when compact enough."""

        residue_ixs = set(getattr(atomgroup.residues, "ix", []))
        if not residue_ixs or len(residue_ixs) > MAX_RESOLVED_OBJECT_REFS:
            return base_object_refs_for_atomgroup(atomgroup)

        refs = [
            self.object_ref_for_residue(residue)
            for residue in atomgroup.universe.residues
            if residue.ix in residue_ixs
        ]
        return refs or base_object_refs_for_atomgroup(atomgroup)

    def annotation_for_residue(self, residue: Any) -> StructuralProfile:
        """Return annotation fields for a residue."""

        profile = self.annotations.residue(residue).to_profile()
        return profile.model_copy(update=self._active_landmarks().profile_update_for_residue(residue.ix))

    def selected_structural_profile(self, atomgroup: Any) -> StructuralProfile | None:
        """Summarize structural annotations for a resolved selection."""

        residues = list(atomgroup.residues)
        if not residues:
            return None
        if len(residues) == 1:
            return self.annotation_for_residue(residues[0])

        secondary_counts: dict[str, int] = {}
        exposure_counts: dict[str, int] = {}
        local_packing_counts: dict[str, int] = {}
        interface_chains: set[str] = set()
        ligand_contacts: set[str] = set()
        landmark_memberships: set[str] = set()
        secondary_elements: set[str] = set()
        ligand_contact_shell_ids: set[str] = set()
        pocket_ids: set[str] = set()
        interchain_contact_interface_ids: set[str] = set()
        biological_interface_ids: set[str] = set()
        definition_ids: set[str] = set()
        reference_keys: set[str] = set()
        for residue in residues:
            annotation = self.annotations.residue(residue)
            landmark_update = self._active_landmarks().profile_update_for_residue(residue.ix)
            annotation_profile = annotation.to_profile()
            if annotation.secondary_structure:
                secondary_counts[annotation.secondary_structure] = (
                    secondary_counts.get(annotation.secondary_structure, 0) + 1
                )
            if annotation.exposure:
                exposure_counts[annotation.exposure] = exposure_counts.get(annotation.exposure, 0) + 1
            if annotation.local_packing:
                local_packing_counts[annotation.local_packing] = (
                    local_packing_counts.get(annotation.local_packing, 0) + 1
                )
            interface_chains.update(annotation.interface_chains)
            ligand_contacts.update(annotation.ligand_contact_ids)
            landmark_memberships.update(landmark_update["landmark_memberships"])
            secondary_element = landmark_update["secondary_structure_element"]
            if secondary_element:
                secondary_elements.add(secondary_element)
            ligand_contact_shell_ids.update(landmark_update["ligand_contact_shell_ids"])
            pocket_ids.update(landmark_update["pocket_ids"])
            interchain_contact_interface_ids.update(
                landmark_update["interchain_contact_interface_ids"]
            )
            biological_interface_ids.update(landmark_update["biological_interface_ids"])
            definition_ids.update(annotation_profile.definition_ids)
            reference_keys.update(annotation_profile.reference_keys)

        return StructuralProfile(
            selected_residue_count=len(residues),
            secondary_structure_counts=dict(sorted(secondary_counts.items())),
            exposure_counts=dict(sorted(exposure_counts.items())),
            local_packing_counts=dict(sorted(local_packing_counts.items())),
            interface_chains=sorted(interface_chains),
            interface_distance_cutoff_A=INTERFACE_DISTANCE_A if interface_chains else None,
            ligand_contacts=sorted(ligand_contacts),
            ligand_contact_cutoff_A=LIGAND_CONTACT_DISTANCE_A if ligand_contacts else None,
            landmark_memberships=sorted(landmark_memberships),
            secondary_structure_element=next(iter(sorted(secondary_elements)), None),
            ligand_contact_shell_ids=sorted(ligand_contact_shell_ids),
            pocket_ids=sorted(pocket_ids),
            interchain_contact_interface_ids=sorted(interchain_contact_interface_ids),
            biological_interface_ids=sorted(biological_interface_ids),
            definition_ids=sorted(definition_ids),
            reference_keys=sorted(reference_keys),
            limitations=list(self.annotations.limitations),
        )

    def nearest_entries(
        self,
        selected: Any,
        object_types: set[str],
        limit: int,
    ) -> list[NearbyObject]:
        """Return nearest residue-like objects for locate-style summaries."""

        selected_ixs = set(selected.residues.ix)
        entries: list[NearbyObject] = []
        for residue in self.universe.residues:
            if residue.ix in selected_ixs:
                continue
            ref = object_ref_for_residue(residue)
            if ref.type not in object_types:
                continue
            distance = min_distance_A(selected, residue.atoms)
            if distance is None:
                continue
            entries.append(
                NearbyObject(
                    object_id=ref.id,
                    min_distance_A=distance,
                    annotations=self.annotation_for_residue(residue),
                )
            )
        entries.sort(key=lambda entry: entry.min_distance_A)
        return entries[:limit]

    def nearby_entries(
        self,
        selected: Any,
        radius_A: float,
        limit: int,
        focus: tuple[ContextFocus, ...] = DEFAULT_CONTEXT_FOCUS,
    ) -> tuple[list[NearbyEntry], bool]:
        """Return ranked nearby residue-like entries for context retrieval."""

        selected_refs = {residue.ix for residue in selected.residues}
        selected_chains = _chains_for_atomgroup(selected)
        selected_heavy = heavy_atomgroup(selected)
        selected_heavy_atoms = list(selected_heavy)
        selected_heavy_positions = selected_heavy.positions
        search_water_bridges = _should_search_water_bridges(selected, focus)
        entries: list[NearbyEntry] = []
        for residue in self.universe.residues:
            if residue.ix in selected_refs:
                continue
            object_type = residue_object_type(residue)
            if object_type == "water":
                continue
            pair = _closest_pair_from_precomputed(
                selected_heavy_atoms,
                selected_heavy_positions,
                residue.atoms,
            )
            if pair is None or pair.distance_A > radius_A:
                continue
            ref = self.object_ref_for_residue(residue)
            relation = self._relation_for_pair(
                selected,
                residue,
                ref.id,
                pair,
                use_backend=bool(set(focus).intersection(INTERACTION_BACKEND_FOCUS)),
            )
            water_relation = (
                self._water_bridge_relation(selected, residue, ref.id, pair)
                if search_water_bridges
                else None
            )
            if water_relation is not None and (
                "water_bridges" in focus or relation.type in {"near", "nonbonded_contact"}
            ):
                relation = water_relation
            residue_chain = chain_for_residue(residue)
            entries.append(
                NearbyEntry(
                    object=ref,
                    relation=relation,
                    distance_A=pair.distance_A,
                    is_interchain=bool(selected_chains and residue_chain not in selected_chains),
                    is_ligand_or_ion=object_type in {"ligand", "ion"},
                    annotation=ref.annotations,
                )
            )

        entries.sort(key=lambda entry: _ranking_key(entry, focus))
        truncated = len(entries) > limit
        return entries[:limit], truncated

    def relation_between(
        self,
        source_atoms: Any,
        target_atoms: Any,
        source_id: str,
        target_id: str,
    ) -> Relation | None:
        """Return closest-pair relation evidence between two selections."""

        pair = closest_heavy_atom_pair(source_atoms, target_atoms)
        if pair is None:
            return None
        return relation_for_atomgroups(source_atoms, target_atoms, source_id, target_id, pair)

    def _relation_for_pair(
        self,
        selected: Any,
        residue: Any,
        target_id: str,
        pair: Any,
        use_backend: bool = False,
    ) -> Relation:
        fallback = relation_for_atomgroups(selected, residue.atoms, "selection", target_id, pair)
        if not use_backend:
            return fallback
        return backend_relation_or_fallback(
            self.interaction_store,
            selected,
            residue.atoms,
            "selection",
            target_id,
            fallback,
        )

    def _water_bridge_relation(
        self,
        selected: Any,
        residue: Any,
        target_id: str,
        direct_pair: Any,
    ) -> Relation | None:
        best: tuple[float, Any, Any, str] | None = None
        for water in self.universe.residues:
            if residue_object_type(water) != "water":
                continue
            source_water_pair = closest_heavy_atom_pair(selected, water.atoms)
            target_water_pair = closest_heavy_atom_pair(residue.atoms, water.atoms)
            if source_water_pair is None or target_water_pair is None:
                continue
            if not _water_bridge_atom_types(source_water_pair, target_water_pair):
                continue
            if (
                not _water_bridge_distance_is_valid(source_water_pair.distance_A)
                or not _water_bridge_distance_is_valid(target_water_pair.distance_A)
            ):
                continue
            score = max(source_water_pair.distance_A, target_water_pair.distance_A)
            water_id = object_ref_for_residue(water).id
            if best is None or score < best[0]:
                best = (score, source_water_pair, target_water_pair, water_id)
        if best is None:
            return None
        _, source_water_pair, target_water_pair, water_id = best
        return water_bridge_relation(
            source_id="selection",
            target_id=target_id,
            direct_pair=direct_pair,
            source_water_pair=source_water_pair,
            target_water_pair=target_water_pair,
            water_id=water_id,
        )


def normalize_context_focus(focus: str | Sequence[str] | None) -> tuple[ContextFocus, ...]:
    """Validate public context focus values."""

    if focus is None:
        return DEFAULT_CONTEXT_FOCUS
    raw_values = [focus] if isinstance(focus, str) else list(focus)
    if not raw_values:
        return DEFAULT_CONTEXT_FOCUS

    allowed = set(CONTEXT_FOCUS_VALUES)
    normalized: list[ContextFocus] = []
    for raw_value in raw_values:
        if not isinstance(raw_value, str):
            raise MetricError(CONTEXT_FOCUS_NOTATION_HELP)
        value = raw_value.strip().lower()
        if value not in allowed:
            raise MetricError(
                f"Unsupported context focus {raw_value!r}. {CONTEXT_FOCUS_NOTATION_HELP}"
            )
        normalized.append(cast(ContextFocus, value))
    return tuple(dict.fromkeys(normalized))


def _non_empty_landmark_updates(updates: dict[str, Any]) -> dict[str, Any]:
    return {key: value for key, value in updates.items() if value}


def _ranking_key(
    entry: NearbyEntry,
    focus: tuple[ContextFocus, ...],
) -> tuple[int, int, int, float, str]:
    focus_priority = 0
    if "interchain_interfaces" in focus and entry.is_interchain:
        focus_priority = -1
    elif "ligand_contact_shell" in focus and entry.is_ligand_or_ion:
        focus_priority = -1
    elif "metal_coordination" in focus and entry.relation.type == "metal_coordination":
        focus_priority = -1
    elif "hydrogen_bonds" in focus and entry.relation.type in {"hydrogen_bond", "polar_contact_candidate"}:
        focus_priority = -1
    elif "salt_bridges" in focus and entry.relation.type == "salt_bridge":
        focus_priority = -1
    elif "hydrophobic_contacts" in focus and entry.relation.type == "hydrophobic_contact":
        focus_priority = -1
    elif "pi_stacking" in focus and entry.relation.type == "pi_stacking":
        focus_priority = -1
    elif "water_bridges" in focus and entry.relation.type == "water_bridge_candidate":
        focus_priority = -1
    elif "steric_clashes" in focus and entry.relation.type == "steric_clash":
        focus_priority = -1
    elif "contacts" in focus and entry.relation.type != "near":
        focus_priority = -1

    relation_priority = _context_relation_priority(entry.relation.type)
    interchain_priority = 0 if entry.is_interchain else 1
    return (focus_priority, relation_priority, interchain_priority, entry.distance_A, entry.object.id)


def _context_relation_priority(relation_type: str) -> int:
    return CONTEXT_RELATION_PRIORITY.get(relation_type, relation_priority_for_type(relation_type))


def _closest_pair_from_precomputed(
    source_atoms: list[Any],
    source_positions: Any,
    target_atoms: Any,
) -> AtomPairDistance | None:
    if not source_atoms:
        return None
    target_heavy = heavy_atomgroup(target_atoms)
    if len(target_heavy) == 0:
        return None

    target_atom_list = list(target_heavy)
    deltas = source_positions[:, None, :] - target_heavy.positions[None, :, :]
    distances = (deltas * deltas).sum(axis=2) ** 0.5
    flat_index = int(distances.argmin())
    source_index, target_index = divmod(flat_index, distances.shape[1])
    return AtomPairDistance(
        distance_A=round(float(distances[source_index, target_index]), 3),
        atom_a=source_atoms[source_index],
        atom_b=target_atom_list[target_index],
    )


def _should_search_water_bridges(selected: Any, focus: tuple[ContextFocus, ...]) -> bool:
    if "water_bridges" in focus:
        return True
    return len(selected.residues) <= 5 and len(selected) <= 200


def _chains_for_atomgroup(atomgroup: Any) -> set[str]:
    return {chain_identifier_for_atom(atom) for atom in atomgroup if chain_identifier_for_atom(atom)}


def _water_bridge_atom_types(source_water_pair: Any, target_water_pair: Any) -> bool:
    return (
        atom_element(source_water_pair.atom_a) in HBOND_DONOR_ACCEPTOR_ELEMENTS
        and atom_element(source_water_pair.atom_b) == "O"
        and atom_element(target_water_pair.atom_a) in HBOND_DONOR_ACCEPTOR_ELEMENTS
        and atom_element(target_water_pair.atom_b) == "O"
    )


def _water_bridge_distance_is_valid(distance_A: float) -> bool:
    return WATER_BRIDGE_MIN_DISTANCE_A <= distance_A <= WATER_BRIDGE_MAX_DISTANCE_A
