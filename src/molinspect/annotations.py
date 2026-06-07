"""Structural annotations computed behind the high-level inspection API."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Any

import numpy as np
from scipy.spatial import cKDTree

from .backends.annotations import (
    SasaProfile as _SasaProfile,
    SecondaryStructureAnnotation as _SecondaryStructure,
    dssp_annotations as _dssp_annotations,
    dssp_key_for_residue as _dssp_key_for_residue,
    freesasa_profiles as _freesasa_profiles,
)
from .objects import object_id_for_residue, object_ref_for_residue, residue_object_type
from .schemas import ExposureSource, LocalPackingSource, StructuralProfile
from .definitions import (
    EXPOSURE_DENSITY_LOWER_PERCENTILE,
    EXPOSURE_DENSITY_UPPER_PERCENTILE,
    EXPOSURE_NEIGHBOR_RADIUS_A,
    INTERFACE_DISTANCE_A,
    LIGAND_CONTACT_DISTANCE_A,
    definition,
)


@dataclass(frozen=True, slots=True)
class ResidueAnnotation:
    """Compact structural profile for one residue-like object."""

    object_id: str
    secondary_structure: str | None = None
    secondary_structure_code: str | None = None
    secondary_structure_source: str | None = None
    exposure: str | None = None
    exposure_source: ExposureSource | None = None
    surface_status: str | None = None
    sasa_A2: float | None = None
    relative_sasa: float | None = None
    local_packing: str | None = None
    local_packing_source: LocalPackingSource | None = None
    local_contact_count: int | None = None
    interface_chains: tuple[str, ...] = ()
    nearest_interchain_distance_A: float | None = None
    ligand_contact_ids: tuple[str, ...] = ()
    limitations: tuple[str, ...] = ()

    @property
    def is_interface(self) -> bool:
        return bool(self.interface_chains)

    def to_profile(self) -> StructuralProfile:
        """Return an LLM-safe annotation profile."""

        return StructuralProfile(
            secondary_structure=self.secondary_structure,
            secondary_structure_code=self.secondary_structure_code,
            secondary_structure_source=self.secondary_structure_source,
            exposure=self.exposure,
            exposure_source=self.exposure_source,
            surface_status=self.surface_status,
            sasa_A2=self.sasa_A2,
            relative_sasa=self.relative_sasa,
            local_packing=self.local_packing,
            local_packing_source=self.local_packing_source,
            local_contact_count=self.local_contact_count,
            local_contact_radius_A=EXPOSURE_NEIGHBOR_RADIUS_A
            if self.local_contact_count is not None
            else None,
            interface_chains=list(self.interface_chains),
            interface_distance_cutoff_A=INTERFACE_DISTANCE_A if self.interface_chains else None,
            nearest_interchain_distance_A=self.nearest_interchain_distance_A,
            ligand_contacts=list(self.ligand_contact_ids),
            ligand_contact_cutoff_A=LIGAND_CONTACT_DISTANCE_A if self.ligand_contact_ids else None,
            definition_ids=_annotation_definition_ids(self),
            reference_keys=_annotation_reference_keys(self),
            limitations=list(self.limitations),
        )

    def to_output(self) -> dict[str, Any]:
        """Return an annotation dictionary suitable for compact object refs."""

        return self.to_profile().model_dump(exclude_none=True, exclude_defaults=True)


@dataclass(slots=True)
class AnnotationStore:
    """Lazily built annotations for the currently loaded topology."""

    by_residue_ix: dict[int, ResidueAnnotation]
    limitations: tuple[str, ...] = ()

    @classmethod
    def build(cls, universe: Any, source_files: tuple[Path, ...]) -> AnnotationStore:
        """Build cheap structural annotations without changing the public API."""

        protein_ixs = _selected_residue_ixs(universe, "protein")
        nucleic_ixs = _selected_residue_ixs(universe, "nucleic")
        classification_sets = (protein_ixs, nucleic_ixs)
        secondary_structure, secondary_limitations = _dssp_annotations(source_files)
        sasa_profiles, sasa_limitations = _freesasa_profiles(source_files)
        local_packing_profiles = _local_packing_profiles(universe, protein_ixs)
        interface_profiles = _interface_profiles(universe, protein_ixs)
        ligand_contacts = _ligand_contact_profiles(universe, classification_sets)

        annotations: dict[int, ResidueAnnotation] = {}
        limitations = [
            definition("interchain_contact_interface").limitations[0],
        ]
        limitations.extend(secondary_limitations)
        limitations.extend(sasa_limitations)

        for residue in universe.residues:
            key = _dssp_key_for_residue(residue)
            secondary = secondary_structure.get(key)
            sasa = sasa_profiles.get(key)
            local_packing = local_packing_profiles.get(residue.ix)
            interface = interface_profiles.get(residue.ix, _InterfaceProfile())
            contacts = ligand_contacts.get(residue.ix, ())
            exposure = sasa.exposure if sasa else None
            exposure_source: ExposureSource | None = "freesasa" if sasa else None
            annotations[residue.ix] = ResidueAnnotation(
                object_id=object_id_for_residue(residue, classification_sets),
                secondary_structure=secondary.name if secondary else None,
                secondary_structure_code=secondary.code if secondary else None,
                secondary_structure_source="mkdssp" if secondary else None,
                exposure=exposure,
                exposure_source=exposure_source,
                surface_status=_surface_status(exposure),
                sasa_A2=sasa.sasa_A2 if sasa else None,
                relative_sasa=sasa.relative_sasa if sasa else None,
                local_packing=local_packing.packing if local_packing else None,
                local_packing_source="local_packing_density" if local_packing else None,
                local_contact_count=local_packing.local_contact_count if local_packing else None,
                interface_chains=interface.chains,
                nearest_interchain_distance_A=interface.nearest_distance_A,
                ligand_contact_ids=contacts,
                limitations=tuple(
                    limitations_for_residue(
                        residue=residue,
                        secondary=secondary,
                        local_packing=local_packing,
                        sasa=sasa,
                        secondary_limitations=secondary_limitations,
                        sasa_limitations=sasa_limitations,
                    )
                ),
            )

        return cls(by_residue_ix=annotations, limitations=tuple(dict.fromkeys(limitations)))

    def residue(self, residue: Any) -> ResidueAnnotation:
        """Return the annotation profile for a residue, or an empty profile."""

        return self.by_residue_ix.get(
            residue.ix,
            ResidueAnnotation(object_id=object_id_for_residue(residue)),
        )

    def object_ref_with_annotations(self, residue: Any) -> Any:
        """Return an object reference with annotation fields attached."""

        annotation = self.residue(residue).to_output()
        ref = object_ref_for_residue(residue)
        return ref.model_copy(update={"annotations": {**ref.annotations, **annotation}})


@dataclass(frozen=True, slots=True)
class _LocalPackingProfile:
    packing: str
    local_contact_count: int


@dataclass(frozen=True, slots=True)
class _InterfaceProfile:
    chains: tuple[str, ...] = ()
    nearest_distance_A: float | None = None


def limitations_for_residue(
    residue: Any,
    secondary: _SecondaryStructure | None,
    local_packing: _LocalPackingProfile | None,
    sasa: _SasaProfile | None,
    secondary_limitations: list[str],
    sasa_limitations: list[str],
) -> list[str]:
    """Return only the annotation limitations relevant to a residue."""

    limitations: list[str] = []
    if residue_object_type(residue) == "residue" and secondary is None:
        limitations.extend(secondary_limitations or ["DSSP secondary structure is unavailable for this residue."])
    if local_packing is not None and sasa is None:
        limitations.extend(sasa_limitations)
        limitations.append(definition("local_packing_density").limitations[0])
    return limitations


def _annotation_definition_ids(annotation: ResidueAnnotation) -> list[str]:
    ids: list[str] = []
    if annotation.secondary_structure:
        ids.append("loop" if annotation.secondary_structure == "loop" else "secondary_structure")
    if annotation.exposure_source == "freesasa":
        ids.append("freesasa_exposure")
    if annotation.local_packing_source == "local_packing_density":
        ids.append("local_packing_density")
    if annotation.interface_chains:
        ids.append("interchain_contact_interface")
    if annotation.ligand_contact_ids:
        ids.append("ligand_contact_shell")
    return list(dict.fromkeys(ids))


def _annotation_reference_keys(annotation: ResidueAnnotation) -> list[str]:
    keys: list[str] = []
    for definition_id in _annotation_definition_ids(annotation):
        keys.extend(definition(definition_id).reference_keys)
    return list(dict.fromkeys(keys))


def _selected_residue_ixs(universe: Any, selection: str) -> set[int]:
    try:
        return set(universe.select_atoms(selection).residues.ix)
    except Exception:
        return set()


def _local_packing_profiles(universe: Any, protein_ixs: set[int]) -> dict[int, _LocalPackingProfile]:
    protein_residues = [residue for residue in universe.residues if residue.ix in protein_ixs]
    if not protein_residues:
        return {}

    centers = np.asarray([_residue_center(residue) for residue in protein_residues], dtype=float)
    deltas = centers[:, None, :] - centers[None, :, :]
    distances = np.sqrt(np.sum(deltas * deltas, axis=2))
    neighbor_counts = np.sum(distances <= EXPOSURE_NEIGHBOR_RADIUS_A, axis=1) - 1
    lower, upper = np.percentile(
        neighbor_counts,
        [EXPOSURE_DENSITY_LOWER_PERCENTILE, EXPOSURE_DENSITY_UPPER_PERCENTILE],
    )

    profiles: dict[int, _LocalPackingProfile] = {}
    for residue, count in zip(protein_residues, neighbor_counts, strict=True):
        count_int = int(count)
        if count <= lower:
            packing = "low_local_packing"
        elif count >= upper:
            packing = "high_local_packing"
        else:
            packing = "medium_local_packing"
        profiles[residue.ix] = _LocalPackingProfile(packing=packing, local_contact_count=count_int)
    return profiles


def _surface_status(exposure: str | None) -> str | None:
    if exposure == "exposed":
        return "surface"
    if exposure == "buried":
        return "core"
    if exposure == "partially_buried":
        return "intermediate"
    return None


def _interface_profiles(universe: Any, protein_ixs: set[int]) -> dict[int, _InterfaceProfile]:
    protein_residues = [residue for residue in universe.residues if residue.ix in protein_ixs]
    positions, residue_ixs, chains, _ = _heavy_atom_table(protein_residues)
    if len(positions) == 0:
        return {}

    chain_contacts: dict[int, set[str]] = {}
    nearest_distances: dict[int, float] = {}
    tree = cKDTree(positions)
    for atom_i, atom_j in tree.query_pairs(r=INTERFACE_DISTANCE_A):
        residue_i = residue_ixs[atom_i]
        residue_j = residue_ixs[atom_j]
        chain_i = chains[atom_i]
        chain_j = chains[atom_j]
        if residue_i == residue_j or chain_i == chain_j:
            continue
        distance = round(float(np.linalg.norm(positions[atom_i] - positions[atom_j])), 3)
        chain_contacts.setdefault(residue_i, set()).add(chain_j)
        chain_contacts.setdefault(residue_j, set()).add(chain_i)
        nearest_distances[residue_i] = min(distance, nearest_distances.get(residue_i, distance))
        nearest_distances[residue_j] = min(distance, nearest_distances.get(residue_j, distance))

    return {
        residue_ix: _InterfaceProfile(
            chains=tuple(sorted(contact_chains)),
            nearest_distance_A=nearest_distances.get(residue_ix),
        )
        for residue_ix, contact_chains in chain_contacts.items()
    }


def _ligand_contact_profiles(
    universe: Any,
    classification_sets: tuple[set[int], set[int]],
) -> dict[int, tuple[str, ...]]:
    protein_residues = [
        residue
        for residue in universe.residues
        if residue_object_type(residue, classification_sets) == "residue"
    ]
    ligand_residues = [
        residue
        for residue in universe.residues
        if residue_object_type(residue, classification_sets) in {"ligand", "ion"}
    ]
    protein_positions, protein_residue_ixs, _, _ = _heavy_atom_table(protein_residues)
    ligand_positions, ligand_residue_ixs, _, ligand_residues_by_ix = _heavy_atom_table(ligand_residues)
    if len(protein_positions) == 0 or len(ligand_positions) == 0:
        return {}

    protein_tree = cKDTree(protein_positions)
    ligand_tree = cKDTree(ligand_positions)
    contacts: dict[int, set[str]] = {}
    for protein_atom_index, ligand_atom_indices in enumerate(
        protein_tree.query_ball_tree(ligand_tree, r=LIGAND_CONTACT_DISTANCE_A)
    ):
        if not ligand_atom_indices:
            continue
        protein_residue_ix = protein_residue_ixs[protein_atom_index]
        for ligand_atom_index in ligand_atom_indices:
            ligand_residue_ix = ligand_residue_ixs[ligand_atom_index]
            ligand_residue = ligand_residues_by_ix[ligand_residue_ix]
            contacts.setdefault(protein_residue_ix, set()).add(
                object_id_for_residue(ligand_residue, classification_sets)
            )
    return {key: tuple(sorted(value)) for key, value in contacts.items()}


def _chain_for_residue(residue: Any) -> str:
    for atom in residue.atoms:
        chain = str(getattr(atom, "chainID", "")).strip()
        if chain:
            return chain
    segment = str(getattr(residue, "segid", "")).strip()
    return segment or "_"


def _residue_center(residue: Any) -> np.ndarray:
    heavy_positions = [
        atom.position for atom in residue.atoms if not str(getattr(atom, "name", "")).upper().startswith("H")
    ]
    if not heavy_positions:
        heavy_positions = [atom.position for atom in residue.atoms]
    return np.asarray(heavy_positions, dtype=float).mean(axis=0)


def _heavy_atom_table(
    residues: list[Any],
) -> tuple[np.ndarray, list[int], list[str], dict[int, Any]]:
    positions: list[np.ndarray] = []
    residue_ixs: list[int] = []
    chains: list[str] = []
    residues_by_ix: dict[int, Any] = {}
    for residue in residues:
        residues_by_ix[residue.ix] = residue
        residue_chain = _chain_for_residue(residue)
        for atom in residue.atoms:
            if str(getattr(atom, "name", "")).upper().startswith("H"):
                continue
            positions.append(np.asarray(atom.position, dtype=float))
            residue_ixs.append(residue.ix)
            chains.append(residue_chain)
    if not positions:
        return np.empty((0, 3), dtype=float), [], [], residues_by_ix
    return np.asarray(positions, dtype=float), residue_ixs, chains, residues_by_ix
