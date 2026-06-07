"""Pydantic schemas for MolInspect outputs."""

from __future__ import annotations

from typing import Any, Literal

from pydantic import BaseModel, ConfigDict, Field


class MolInspectModel(BaseModel):
    model_config = ConfigDict(extra="forbid")


class EvidenceItem(MolInspectModel):
    type: str
    metric: str | None = None
    value: Any | None = None
    unit: str | None = None
    frame: int | None = None
    source: str | None = None


ObjectType = Literal[
    "atom",
    "residue",
    "ligand",
    "ion",
    "water",
    "chain",
    "selection_region",
    "secondary_structure",
    "loop",
    "ligand_contact_shell",
    "pocket",
    "interchain_contact_interface",
    "biological_interface",
]
ContextFocus = Literal[
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
]
ContextBudget = Literal["small", "medium", "large"]
ContextScale = Literal[
    "chemical_contacts",
    "residue_environment",
    "ligand_binding_site",
    "metal_coordination",
    "protein_interface",
    "broad_environment",
]
ExposureSource = Literal["freesasa"]
LocalPackingSource = Literal["local_packing_density"]


class StructuralProfile(MolInspectModel):
    selected_residue_count: int | None = None
    secondary_structure: str | None = None
    secondary_structure_code: str | None = None
    secondary_structure_source: str | None = None
    secondary_structure_counts: dict[str, int] = Field(default_factory=dict)
    exposure: str | None = None
    exposure_source: ExposureSource | None = None
    surface_status: str | None = None
    sasa_A2: float | None = None
    relative_sasa: float | None = None
    local_packing: str | None = None
    local_packing_source: LocalPackingSource | None = None
    local_contact_count: int | None = None
    local_contact_radius_A: float | None = None
    exposure_counts: dict[str, int] = Field(default_factory=dict)
    local_packing_counts: dict[str, int] = Field(default_factory=dict)
    interface_chains: list[str] = Field(default_factory=list)
    interface_distance_cutoff_A: float | None = None
    nearest_interchain_distance_A: float | None = None
    ligand_contacts: list[str] = Field(default_factory=list)
    ligand_contact_cutoff_A: float | None = None
    landmark_memberships: list[str] = Field(default_factory=list)
    secondary_structure_element: str | None = None
    ligand_contact_shell_ids: list[str] = Field(default_factory=list)
    pocket_ids: list[str] = Field(default_factory=list)
    interchain_contact_interface_ids: list[str] = Field(default_factory=list)
    biological_interface_ids: list[str] = Field(default_factory=list)
    definition_ids: list[str] = Field(default_factory=list)
    reference_keys: list[str] = Field(default_factory=list)
    limitations: list[str] = Field(default_factory=list)


class NearbyObject(MolInspectModel):
    object_id: str
    min_distance_A: float
    annotations: StructuralProfile | None = None


class ObjectRef(MolInspectModel):
    id: str
    type: ObjectType
    name: str | None = None
    chain: str | None = None
    segment: str | None = None
    resid: int | str | None = None
    icode: str | None = None
    resname: str | None = None
    atom_name: str | None = None
    altloc: str | None = None
    atom_index: int | None = None
    atom_count: int | None = None
    annotations: dict[str, Any] = Field(default_factory=dict)


class Relation(MolInspectModel):
    source: str
    target: str
    type: str
    category: str | None = None
    confidence: str | None = None
    backend: str | None = None
    definition_id: str | None = None
    definition_source: str | None = None
    reference_keys: list[str] = Field(default_factory=list)
    min_distance_A: float | None = None
    cutoff_A: float | None = None
    angle_deg: float | None = None
    source_atom: str | None = None
    target_atom: str | None = None
    occupancy: float | None = None
    frames: list[int] | None = None
    evidence: list[EvidenceItem] = Field(default_factory=list)
    limitations: list[str] = Field(default_factory=list)


class LoadResult(MolInspectModel):
    target_id: str
    name: str | None = None
    n_atoms: int
    n_residues: int
    n_chains: int | None = None
    n_segments: int | None = None
    n_frames: int
    mode: Literal["static", "trajectory"]
    source_files: list[str] = Field(default_factory=list)
    available_annotations: list[str] = Field(default_factory=list)
    limitations: list[str] = Field(default_factory=list)


class ObjectsResult(MolInspectModel):
    objects: list[ObjectRef]
    count: int
    returned: int
    truncated: bool = False
    limitations: list[str] = Field(default_factory=list)


class SelectionResult(MolInspectModel):
    selection: str
    expression: str
    frame: int
    n_atoms: int
    resolved_objects: list[str]
    limitations: list[str] = Field(default_factory=list)


class LocationInfo(MolInspectModel):
    center_of_geometry_A: list[float]
    selected_atom_count: int
    selected_object_count: int
    chain: str | list[str] | None = None
    segment: str | list[str] | None = None
    objects: list[str] = Field(default_factory=list)
    structural_profile: StructuralProfile | None = None
    secondary_structure: str | None = None
    exposure_status: str | None = None
    exposure_source: ExposureSource | None = None
    surface_status: str | None = None
    sasa_A2: float | None = None
    relative_sasa: float | None = None
    local_packing: str | None = None
    local_packing_source: LocalPackingSource | None = None
    local_contact_count: int | None = None
    local_contact_radius_A: float | None = None
    interface_chains: list[str] = Field(default_factory=list)
    nearest_interchain_distance_A: float | None = None
    ligand_contacts: list[str] = Field(default_factory=list)
    landmark_memberships: list[str] = Field(default_factory=list)
    secondary_structure_element: str | None = None
    ligand_contact_shell_ids: list[str] = Field(default_factory=list)
    pocket_ids: list[str] = Field(default_factory=list)
    interchain_contact_interface_ids: list[str] = Field(default_factory=list)
    biological_interface_ids: list[str] = Field(default_factory=list)
    nearest_protein_residues: list[NearbyObject] = Field(default_factory=list)
    near_ligands: list[NearbyObject] = Field(default_factory=list)
    distance_to_chain_centroid_A: float | None = None
    plain_language: str


class LocateResult(MolInspectModel):
    selection: str
    resolved_objects: list[str]
    frame: int
    location: LocationInfo
    evidence: list[EvidenceItem] = Field(default_factory=list)
    limitations: list[str] = Field(default_factory=list)


class ContextResult(MolInspectModel):
    selection: str
    frame: int
    scale: ContextScale | None = None
    radius_A: float
    budget: ContextBudget = "small"
    focus: list[ContextFocus] = Field(default_factory=list)
    objects: list[ObjectRef]
    relations: list[Relation]
    summary: str
    evidence: list[EvidenceItem] = Field(default_factory=list)
    limitations: list[str] = Field(default_factory=list)


class ContextScaleInfo(MolInspectModel):
    id: ContextScale
    radius_A: float
    budget: ContextBudget
    focus: list[ContextFocus]
    description: str


class ContextScalesResult(MolInspectModel):
    scales: list[ContextScaleInfo]
    count: int
    limitations: list[str] = Field(default_factory=list)


class TimelineMetricInfo(MolInspectModel):
    metric: str
    required_arguments: list[str]
    optional_arguments: list[str] = Field(default_factory=list)
    best_for: str
    description: str
    summary_keys: list[str]
    sampled_value_keys: list[str]
    event_types: list[str]
    example: str
    limitations: list[str] = Field(default_factory=list)


class TimelineMetricsResult(MolInspectModel):
    metrics: list[TimelineMetricInfo]
    count: int
    limitations: list[str] = Field(default_factory=list)


class TimelineResult(MolInspectModel):
    metric: str
    selection: str | None = None
    selection1: str | None = None
    selection2: str | None = None
    frames_analyzed: int
    summary_text: str = ""
    summary: dict[str, Any]
    events: list[dict[str, Any]] = Field(default_factory=list)
    sampled_values: list[dict[str, Any]] = Field(default_factory=list)
    limitations: list[str] = Field(default_factory=list)


class CompareResult(MolInspectModel):
    selection: str
    frame_a: int
    frame_b: int
    main_changes: list[dict[str, Any]]
    summary: str
    evidence: list[EvidenceItem] = Field(default_factory=list)
    limitations: list[str] = Field(default_factory=list)
