"""Canonical structural definitions used by MolInspect internals.

This module owns MolInspect's structural vocabulary. Runtime code should import
definition-backed constants from here instead of redefining semantic cutoffs,
labels, source names, or limitation text in local modules.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from types import MappingProxyType
from typing import Any, Literal, Mapping

DefinitionKind = Literal["object", "annotation", "relation", "metric", "temporal_event"]
DefinitionSource = Literal[
    "topology",
    "backend",
    "literature_informed_heuristic",
    "molinspect_policy",
]

MOLINSPECT_HEURISTIC_BACKEND = "molinspect_heuristic"


@dataclass(frozen=True, slots=True)
class Reference:
    """Short reference metadata for definitions.

    The key is what output models and docs expose. The URL is kept here so tests
    and docs can verify that definition references are not invented ad hoc.
    """

    key: str
    label: str
    url: str
    note: str


@dataclass(frozen=True, slots=True)
class StructuralDefinition:
    """One auditable structural term used by MolInspect."""

    id: str
    label: str
    kind: DefinitionKind
    source: DefinitionSource
    description: str
    parameters: Mapping[str, Any] = field(default_factory=dict)
    reference_keys: tuple[str, ...] = ()
    limitations: tuple[str, ...] = ()
    backend: str | None = None
    category: str | None = None
    default_confidence: str | None = None
    priority: int | None = None

    def __post_init__(self) -> None:
        object.__setattr__(self, "parameters", MappingProxyType(dict(self.parameters)))


REFERENCE_CATALOG: dict[str, Reference] = {
    "mdanalysis_universe": Reference(
        key="mdanalysis_universe",
        label="MDAnalysis Universe",
        url="https://docs.mdanalysis.org/stable/documentation_pages/core/universe.html",
        note="MolInspect's required topology/trajectory backend.",
    ),
    "gemmi_docs": Reference(
        key="gemmi_docs",
        label="Gemmi documentation",
        url="https://gemmi.readthedocs.io/en/stable/analysis.html",
        note="Optional static structure parsing and mmCIF/PDB support.",
    ),
    "dssp_docs": Reference(
        key="dssp_docs",
        label="DSSP secondary structure code",
        url="https://swift.cmbi.umcn.nl/gv/dssp/DSSP_2.html",
        note="DSSP code mapping for secondary-structure annotations.",
    ),
    "freesasa_docs": Reference(
        key="freesasa_docs",
        label="FreeSASA Python module",
        url="https://freesasa.github.io/python/classes.html",
        note="Optional exact solvent-accessible surface area backend.",
    ),
    "plip_config": Reference(
        key="plip_config",
        label="PLIP detection thresholds",
        url="https://raw.githubusercontent.com/pharmai/plip/master/plip/basic/config.py",
        note="Default PLIP geometric thresholds for non-covalent interaction detection.",
    ),
    "plip_docs": Reference(
        key="plip_docs",
        label="PLIP documentation",
        url="https://github.com/pharmai/plip/blob/master/DOCUMENTATION.md",
        note="Rule-based non-covalent interaction detection concepts.",
    ),
    "arpeggio_paper": Reference(
        key="arpeggio_paper",
        label="PDBe-Arpeggio paper",
        url="https://pmc.ncbi.nlm.nih.gov/articles/PMC5282402/",
        note="Backend macromolecular interatomic interaction classes.",
    ),
    "p2rank": Reference(
        key="p2rank",
        label="P2Rank",
        url="https://github.com/rdk/p2rank",
        note="Preferred optional ligand-binding site prediction backend.",
    ),
    "fpocket": Reference(
        key="fpocket",
        label="fpocket",
        url="https://github.com/Discngine/fpocket",
        note="Optional geometric pocket fallback when P2Rank is unavailable.",
    ),
    "pisa_docs": Reference(
        key="pisa_docs",
        label="PISA interface analysis",
        url="https://cloud.ccp4.ac.uk/manuals/html-taskref/doc.task.PISA.html",
        note="Reference backend for solvent-surface macromolecular interface analysis.",
    ),
    "mdanalysis_rms": Reference(
        key="mdanalysis_rms",
        label="MDAnalysis RMSD/RMSF documentation",
        url="https://docs.mdanalysis.org/2.4.1/documentation_pages/analysis/rms.html",
        note="RMSD/RMSF definitions and alignment caveats.",
    ),
    "mdanalysis_hbonds": Reference(
        key="mdanalysis_hbonds",
        label="MDAnalysis HydrogenBondAnalysis",
        url="https://docs.mdanalysis.org/stable/documentation_pages/analysis/hydrogenbonds.html",
        note="Trajectory hydrogen-bond detection with donor, hydrogen, acceptor, distance, and angle evidence.",
    ),
    "mdanalysis_distances": Reference(
        key="mdanalysis_distances",
        label="MDAnalysis distance analysis",
        url="https://docs.mdanalysis.org/stable/documentation_pages/analysis/distances.html",
        note="Distance/contact calculations with optional periodic-boundary handling.",
    ),
    "mdanalysis_guesser": Reference(
        key="mdanalysis_guesser",
        label="MDAnalysis topology guessers",
        url="https://docs.mdanalysis.org/stable/documentation_pages/guesser_modules/default_guesser.html",
        note="Distance-based bond and aromaticity guessing conventions.",
    ),
    "molprobity_clash": Reference(
        key="molprobity_clash",
        label="MolProbity clash overlap convention",
        url="https://pmc.ncbi.nlm.nih.gov/articles/PMC2803126/",
        note="Steric clash convention based on non-bonded van der Waals overlap.",
    ),
    "bondi_vdw_radii": Reference(
        key="bondi_vdw_radii",
        label="Bondi van der Waals radii",
        url="https://doi.org/10.1021/j100785a001",
        note="Common van der Waals radius source for steric overlap estimates.",
    ),
    "pdb_ccd": Reference(
        key="pdb_ccd",
        label="PDB Chemical Component Dictionary",
        url="https://www.wwpdb.org/data/ccd",
        note="Preferred future source for residue/ligand/ion identity and covalent component data.",
    ),
}


SECONDARY_STRUCTURE_NAMES = {
    "H": "alpha_helix",
    "B": "beta_bridge",
    "E": "beta_strand",
    "G": "three_ten_helix",
    "I": "pi_helix",
    "T": "turn",
    "S": "bend",
    " ": "loop",
}

WATER_RESNAMES = frozenset({"HOH", "WAT", "H2O", "TIP3", "TIP3P", "SOL", "SOLV", "SPC"})
ION_RESNAMES = frozenset(
    {
        "AL",
        "BA",
        "BR",
        "CA",
        "CD",
        "CL",
        "CO",
        "CS",
        "CU",
        "FE",
        "HG",
        "IOD",
        "K",
        "LI",
        "MG",
        "MN",
        "NA",
        "NI",
        "RB",
        "SR",
        "ZN",
    }
)

METAL_ELEMENTS = {
    "AG",
    "AL",
    "AU",
    "BA",
    "CA",
    "CD",
    "CE",
    "CO",
    "CR",
    "CS",
    "CU",
    "EU",
    "FE",
    "GA",
    "GD",
    "HG",
    "IN",
    "IR",
    "K",
    "LA",
    "LI",
    "LU",
    "MG",
    "MN",
    "NA",
    "NI",
    "OS",
    "PB",
    "PD",
    "PT",
    "RB",
    "RH",
    "RU",
    "SB",
    "SM",
    "SR",
    "TB",
    "TL",
    "W",
    "YB",
    "ZN",
}
COORDINATING_ELEMENTS = {"N", "O", "S", "P", "CL", "BR", "I"}
HBOND_DONOR_ACCEPTOR_ELEMENTS = {"N", "O", "S"}
HYDROPHOBIC_RESNAMES = {
    "ALA",
    "VAL",
    "LEU",
    "ILE",
    "MET",
    "PHE",
    "TRP",
    "PRO",
    "TYR",
}
AROMATIC_RESNAMES = {"PHE", "TYR", "TRP", "HIS"}
AROMATIC_RING_ATOMS = {
    "PHE": ("CG", "CD1", "CD2", "CE1", "CE2", "CZ"),
    "TYR": ("CG", "CD1", "CD2", "CE1", "CE2", "CZ"),
    "TRP": ("CD2", "CE2", "CE3", "CZ2", "CZ3", "CH2"),
    "HIS": ("CG", "ND1", "CD2", "CE1", "NE2"),
}
POSITIVE_ATOMS = {
    "ARG": {"NH1", "NH2", "NE", "CZ"},
    "LYS": {"NZ"},
    "HIS": {"ND1", "NE2"},
}
NEGATIVE_ATOMS = {
    "ASP": {"OD1", "OD2"},
    "GLU": {"OE1", "OE2"},
}


_DEFINITIONS = (
    StructuralDefinition(
        id="atom",
        label="atom",
        kind="object",
        source="topology",
        backend="MDAnalysis",
        description="Atom identity and coordinates as loaded from the molecular topology.",
        reference_keys=("mdanalysis_universe",),
    ),
    StructuralDefinition(
        id="residue",
        label="residue",
        kind="object",
        source="topology",
        backend="MDAnalysis",
        description="Residue-like topology object classified as protein or nucleic by backend selection.",
        reference_keys=("mdanalysis_universe",),
    ),
    StructuralDefinition(
        id="chain",
        label="chain",
        kind="object",
        source="topology",
        backend="MDAnalysis",
        description="PDB chain ID when present; otherwise segment ID, then '_'.",
        reference_keys=("mdanalysis_universe",),
    ),
    StructuralDefinition(
        id="water",
        label="water",
        kind="object",
        source="molinspect_policy",
        description="Residue classified as water by a curated residue-name set.",
        parameters={"resnames": tuple(sorted(WATER_RESNAMES))},
        reference_keys=("pdb_ccd",),
        limitations=("Water identity uses residue-name policy when chemical component metadata is unavailable.",),
    ),
    StructuralDefinition(
        id="ion",
        label="ion",
        kind="object",
        source="molinspect_policy",
        description="Residue classified as ion by residue name and small atom count.",
        parameters={"resnames": tuple(sorted(ION_RESNAMES)), "max_atom_count": 4},
        reference_keys=("pdb_ccd",),
        limitations=("Ion identity uses residue-name policy when chemical component metadata is unavailable.",),
    ),
    StructuralDefinition(
        id="ligand",
        label="ligand",
        kind="object",
        source="molinspect_policy",
        description="Residue that is not protein, nucleic, water, or ion.",
        parameters={"excludes": ("protein", "nucleic", "water", "ion")},
        reference_keys=("mdanalysis_universe", "pdb_ccd"),
        limitations=("Ligand identity is an exclusion rule until chemical component metadata is integrated.",),
    ),
    StructuralDefinition(
        id="selection_region",
        label="selection_region",
        kind="object",
        source="molinspect_policy",
        description="Synthetic compact object for broad residue selections.",
        parameters={"max_residue_objects": 25},
        limitations=("Selection-region IDs are compact handles, not biological domains.",),
    ),
    StructuralDefinition(
        id="secondary_structure",
        label="secondary_structure",
        kind="object",
        source="backend",
        backend="mkdssp",
        description="Consecutive non-loop DSSP secondary-structure run.",
        parameters={"dssp_code_map": SECONDARY_STRUCTURE_NAMES},
        reference_keys=("dssp_docs",),
        limitations=("Secondary-structure objects are emitted only when DSSP annotations are available.",),
    ),
    StructuralDefinition(
        id="loop",
        label="loop",
        kind="object",
        source="backend",
        backend="mkdssp",
        description="Consecutive DSSP loop/irregular run.",
        parameters={"dssp_code": " "},
        reference_keys=("dssp_docs",),
        limitations=("Loop objects are emitted only when DSSP annotations are available.",),
    ),
    StructuralDefinition(
        id="local_packing_density",
        label="local_packing_density",
        kind="annotation",
        source="molinspect_policy",
        description="Residue-center local-packing tier used when exact solvent exposure is unavailable.",
        parameters={"neighbor_radius_A": 10.0, "lower_percentile": 33, "upper_percentile": 66},
        limitations=(
            "Local packing density is not solvent exposure and must not be interpreted as FreeSASA/RSA.",
        ),
    ),
    StructuralDefinition(
        id="freesasa_exposure",
        label="freesasa",
        kind="annotation",
        source="backend",
        backend="FreeSASA",
        description="Residue solvent-accessible surface area from FreeSASA residue areas.",
        parameters={"buried_relative_sasa_max": 0.05, "exposed_relative_sasa_min": 0.25},
        reference_keys=("freesasa_docs",),
        limitations=(
            "Buried/exposed labels use configurable RSA thresholds; inspect relative_sasa for exact values.",
        ),
    ),
    StructuralDefinition(
        id="interchain_contact_interface",
        label="interchain_contact_interface",
        kind="object",
        source="literature_informed_heuristic",
        description="Chain-pair landmark from inter-chain heavy-atom contacts.",
        parameters={"contact_cutoff_A": 5.0, "method": "inter_chain_heavy_atom_contacts"},
        reference_keys=("pisa_docs",),
        limitations=(
            "This is a contact interface, not a PISA biological assembly/interface-energy call.",
        ),
    ),
    StructuralDefinition(
        id="biological_interface",
        label="biological_interface",
        kind="object",
        source="backend",
        backend="PISA",
        description="Macromolecular interface reported by PISA from buried surface and interface energetics.",
        parameters={"method": "pisa_surface_and_assembly_interface_analysis"},
        reference_keys=("pisa_docs",),
        limitations=(
            "Biological-interface objects are emitted only when local PISA runs successfully; "
            "PISA reports should be interpreted with crystallographic assembly context.",
        ),
    ),
    StructuralDefinition(
        id="ligand_contact_shell",
        label="ligand_contact_shell",
        kind="object",
        source="literature_informed_heuristic",
        description="Ligand or ion plus protein residues in its direct heavy-atom contact shell.",
        parameters={"contact_cutoff_A": 4.0, "method": "ligand_contact_shell"},
        reference_keys=("plip_docs", "plip_config"),
        limitations=("This is a direct contact shell, not a full predicted binding site.",),
    ),
    StructuralDefinition(
        id="pocket",
        label="pocket",
        kind="object",
        source="backend",
        backend="P2Rank_or_fpocket",
        description="Geometric or ligandability-scored pocket from a dedicated pocket backend.",
        parameters={"method": "backend_detected_pocket"},
        reference_keys=("p2rank", "fpocket"),
        limitations=("Pocket objects are emitted only when a pocket backend succeeds.",),
    ),
    StructuralDefinition(
        id="topology_bond",
        label="topology_bond",
        kind="relation",
        source="topology",
        backend="MDAnalysis",
        description="Covalent bond read from topology or explicit connectivity records.",
        reference_keys=("mdanalysis_universe",),
        category="covalent",
        default_confidence="topology",
        priority=0,
    ),
    StructuralDefinition(
        id="inferred_covalent_bond",
        label="inferred_covalent_bond",
        kind="relation",
        source="literature_informed_heuristic",
        backend=MOLINSPECT_HEURISTIC_BACKEND,
        description="Covalent bond inferred from local residue/polymer geometry when topology lacks bonds.",
        parameters={
            "vdw_fudge_factor": 0.55,
            "minimum_distance_A": 0.1,
            "peptide_c_n_max_distance_A": 1.7,
        },
        reference_keys=("mdanalysis_guesser",),
        category="covalent",
        default_confidence="inferred",
        priority=1,
        limitations=("Inferred covalent bonds are local connectivity estimates, not force-field topology.",),
    ),
    StructuralDefinition(
        id="steric_clash",
        label="steric_clash",
        kind="relation",
        source="literature_informed_heuristic",
        backend=MOLINSPECT_HEURISTIC_BACKEND,
        description="Non-bonded heavy atoms whose van der Waals spheres overlap beyond the clash threshold.",
        parameters={"min_vdw_overlap_A": 0.4},
        reference_keys=("molprobity_clash", "bondi_vdw_radii"),
        category="steric",
        default_confidence="geometry",
        priority=2,
        limitations=("Steric clashes exclude topology and inferred covalent bonds.",),
    ),
    StructuralDefinition(
        id="metal_coordination",
        label="metal_coordination",
        kind="relation",
        source="literature_informed_heuristic",
        backend=MOLINSPECT_HEURISTIC_BACKEND,
        description="Metal/coordinating atom pair within the coordination distance cutoff.",
        parameters={
            "max_distance_A": 3.0,
            "metal_elements": tuple(sorted(METAL_ELEMENTS)),
            "coordinating_elements": tuple(sorted(COORDINATING_ELEMENTS)),
        },
        reference_keys=("plip_docs", "plip_config", "arpeggio_paper"),
        category="coordination",
        default_confidence="geometry",
        priority=3,
        limitations=("Metal coordination uses distance and element rules, not full valence geometry.",),
    ),
    StructuralDefinition(
        id="salt_bridge",
        label="salt_bridge",
        kind="relation",
        source="literature_informed_heuristic",
        backend=MOLINSPECT_HEURISTIC_BACKEND,
        description="Oppositely charged residue atom names within the salt-bridge cutoff.",
        parameters={
            "max_distance_A": 5.5,
            "positive_atoms": POSITIVE_ATOMS,
            "negative_atoms": NEGATIVE_ATOMS,
        },
        reference_keys=("plip_docs", "plip_config", "arpeggio_paper"),
        category="electrostatic",
        default_confidence="geometry",
        priority=4,
        limitations=(
            "Salt bridges use residue atom-name charge proxies, not full protonation/charge-center inference.",
        ),
    ),
    StructuralDefinition(
        id="hydrogen_bond",
        label="hydrogen_bond",
        kind="relation",
        source="literature_informed_heuristic",
        backend=MOLINSPECT_HEURISTIC_BACKEND,
        description="Donor/acceptor atom pair with explicit hydrogen-angle validation.",
        parameters={
            "donor_acceptor_elements": tuple(sorted(HBOND_DONOR_ACCEPTOR_ELEMENTS)),
            "max_distance_A": 4.1,
            "min_angle_deg": 100.0,
            "donor_hydrogen_max_distance_A": 1.25,
        },
        reference_keys=("plip_docs", "plip_config", "arpeggio_paper"),
        category="polar",
        default_confidence="geometry",
        priority=5,
    ),
    StructuralDefinition(
        id="polar_contact_candidate",
        label="polar_contact_candidate",
        kind="relation",
        source="literature_informed_heuristic",
        backend=MOLINSPECT_HEURISTIC_BACKEND,
        description="Polar atom pair within hydrogen-bond distance without donor/acceptor-angle validation.",
        parameters={
            "donor_acceptor_elements": tuple(sorted(HBOND_DONOR_ACCEPTOR_ELEMENTS)),
            "max_distance_A": 4.1,
        },
        reference_keys=("plip_docs", "plip_config", "arpeggio_paper"),
        category="polar",
        default_confidence="candidate",
        priority=10,
        limitations=("Polar contact candidates lack explicit donor/acceptor and angle validation.",),
    ),
    StructuralDefinition(
        id="pi_stacking",
        label="pi_stacking",
        kind="relation",
        source="literature_informed_heuristic",
        backend=MOLINSPECT_HEURISTIC_BACKEND,
        description="Aromatic ring pair passing center-distance and ring-plane angle checks.",
        parameters={
            "max_center_distance_A": 5.5,
            "max_angle_deviation_deg": 30.0,
            "max_parallel_offset_A": 2.0,
            "aromatic_resnames": tuple(sorted(AROMATIC_RESNAMES)),
            "ring_atoms": AROMATIC_RING_ATOMS,
        },
        reference_keys=("plip_docs", "plip_config", "arpeggio_paper"),
        category="aromatic",
        default_confidence="geometry",
        priority=6,
        limitations=("Protein-residue ring templates are used; ligand aromaticity needs RDKit/PLIP backend.",),
    ),
    StructuralDefinition(
        id="hydrophobic_contact",
        label="hydrophobic_contact",
        kind="relation",
        source="literature_informed_heuristic",
        backend=MOLINSPECT_HEURISTIC_BACKEND,
        description="Hydrophobic side-chain carbon pair within the contact cutoff.",
        parameters={
            "max_distance_A": 4.0,
            "hydrophobic_resnames": tuple(sorted(HYDROPHOBIC_RESNAMES)),
        },
        reference_keys=("plip_docs", "plip_config", "arpeggio_paper"),
        category="hydrophobic",
        default_confidence="geometry",
        priority=7,
        limitations=("Hydrophobic contacts use residue/atom-name proxies, not full ligand atom typing.",),
    ),
    StructuralDefinition(
        id="water_bridge_candidate",
        label="water_bridge_candidate",
        kind="relation",
        source="literature_informed_heuristic",
        backend=MOLINSPECT_HEURISTIC_BACKEND,
        description="One water oxygen bridging two polar atoms by PLIP-like distance bounds.",
        parameters={
            "polar_elements": tuple(sorted(HBOND_DONOR_ACCEPTOR_ELEMENTS)),
            "water_element": "O",
            "min_leg_distance_A": 2.5,
            "max_leg_distance_A": 4.1,
        },
        reference_keys=("plip_docs", "plip_config"),
        category="water_mediated",
        default_confidence="candidate",
        priority=8,
        limitations=(
            "Water bridge candidate lacks PLIP omega/theta angle validation unless hydrogens are available.",
        ),
    ),
    StructuralDefinition(
        id="nonbonded_contact",
        label="nonbonded_contact",
        kind="relation",
        source="literature_informed_heuristic",
        backend=MOLINSPECT_HEURISTIC_BACKEND,
        description="Closest non-covalent heavy-atom distance at or below the generic contact cutoff.",
        parameters={"max_distance_A": 4.0},
        reference_keys=("plip_docs", "plip_config", "arpeggio_paper"),
        category="generic_contact",
        default_confidence="geometry",
        priority=11,
        limitations=("Generic contacts are distance evidence, not chemically typed interactions.",),
    ),
    StructuralDefinition(
        id="near",
        label="near",
        kind="relation",
        source="molinspect_policy",
        backend=MOLINSPECT_HEURISTIC_BACKEND,
        description="Within the requested context radius but outside the generic contact cutoff.",
        parameters={"contact_cutoff_A": 4.0},
        category="proximity",
        default_confidence="geometry",
        priority=12,
    ),
    StructuralDefinition(
        id="centroid_distance",
        label="centroid_distance",
        kind="metric",
        source="molinspect_policy",
        description=(
            "Distance between centers of geometry of two selections over frames, useful for "
            "domain, loop, partner, or site opening/closing motions."
        ),
        parameters={"center": "center_of_geometry"},
        limitations=(
            "Centroid distance is not a contact distance and is not mass-weighted.",
            "Interpretation depends on the selected structural regions.",
        ),
    ),
    StructuralDefinition(
        id="rmsd",
        label="rmsd",
        kind="metric",
        source="backend",
        description="Kabsch-aligned RMSD over the selected atoms against the first selected frame.",
        parameters={"alignment": "kabsch_over_selection"},
        reference_keys=("mdanalysis_rms",),
        limitations=("RMSD requires consistent atom correspondence and whole/aligned molecules.",),
        backend="MDAnalysis.analysis.rms",
    ),
    StructuralDefinition(
        id="rmsf",
        label="rmsf",
        kind="metric",
        source="backend",
        description="Kabsch-aligned per-atom RMSF summarized by selected residue/object.",
        parameters={"alignment": "kabsch_over_selection"},
        reference_keys=("mdanalysis_rms",),
        limitations=("RMSF requires consistent atom correspondence and whole/aligned molecules.",),
        backend="MDAnalysis.analysis.rms",
    ),
    StructuralDefinition(
        id="hydrogen_bond_occupancy",
        label="hydrogen_bond_occupancy",
        kind="metric",
        source="backend",
        description=(
            "Frame occupancy of explicit-hydrogen donor-acceptor hydrogen bonds between two "
            "selected atom groups."
        ),
        parameters={
            "donor_acceptor_distance_cutoff_A": 3.0,
            "donor_hydrogen_distance_cutoff_A": 1.2,
            "donor_hydrogen_acceptor_angle_min_deg": 150.0,
        },
        reference_keys=("mdanalysis_hbonds",),
        limitations=(
            "Requires explicit hydrogens or donor-hydrogen topology that MDAnalysis can resolve.",
            "Occupancy is selection-dependent and does not estimate energetic stability.",
        ),
        backend="MDAnalysis.analysis.hydrogenbonds.HydrogenBondAnalysis",
    ),
    StructuralDefinition(
        id="interaction_persistence",
        label="interaction_persistence",
        kind="metric",
        source="molinspect_policy",
        description=(
            "Combined contact, relation, and explicit-hydrogen hydrogen-bond persistence "
            "summary between two selected atom groups."
        ),
        parameters={"contact_cutoff_A": 4.0},
        reference_keys=("mdanalysis_distances", "mdanalysis_hbonds"),
        limitations=(
            "Non-H-bond relation labels still use MolInspect relation heuristics unless a backend relation is available.",
            "Persistence is selection-dependent and does not estimate binding free energy.",
        ),
    ),
    StructuralDefinition(
        id="ligand_stability",
        label="ligand_stability",
        kind="metric",
        source="molinspect_policy",
        description=(
            "Ligand/site stability summary from site-aligned ligand RMSD, contact occupancy, "
            "and explicit-hydrogen hydrogen-bond occupancy."
        ),
        parameters={
            "stable_min_contact_occupancy": 0.8,
            "stable_max_site_aligned_ligand_rmsd_A": 2.0,
            "unstable_max_contact_occupancy": 0.5,
            "unstable_min_site_aligned_ligand_rmsd_A": 5.0,
        },
        reference_keys=("mdanalysis_rms", "mdanalysis_distances", "mdanalysis_hbonds"),
        limitations=(
            "The stability class is a MolInspect summary rule, not a thermodynamic binding-stability estimate.",
            "Site-aligned ligand RMSD requires stable atom correspondence in ligand and site selections.",
        ),
    ),
    StructuralDefinition(
        id="conformational_states",
        label="conformational_states",
        kind="metric",
        source="molinspect_policy",
        description=(
            "Simple frame-state summary from aligned RMSD and selection-spread features."
        ),
        parameters={
            "min_rmsd_range_A_for_two_states": 1.0,
            "max_state_count": 2,
        },
        reference_keys=("mdanalysis_rms",),
        limitations=(
            "State summaries are deterministic feature bins, not kinetic clustering or Markov-state models.",
            "Representative frames depend on the selected atoms and sampled frame range.",
        ),
    ),
    StructuralDefinition(
        id="selection_spread",
        label="selection_spread",
        kind="metric",
        source="molinspect_policy",
        description="Radius-of-gyration spread of selected atoms across frames.",
        parameters={"spread_proxy": "selection_radius_of_gyration"},
        limitations=("Selection spread is not pocket volume, pocket opening, or ligand accessibility.",),
    ),
    StructuralDefinition(
        id="context_radius",
        label="context_radius",
        kind="metric",
        source="molinspect_policy",
        description="Default local context retrieval radius.",
        parameters={"default_radius_A": 8.0},
    ),
    StructuralDefinition(
        id="context_budget",
        label="context_budget",
        kind="metric",
        source="molinspect_policy",
        description="Named caps for nearby objects returned by context().",
        parameters={"small": 25, "medium": 75, "large": 200},
    ),
    StructuralDefinition(
        id="context_scale",
        label="context_scale",
        kind="metric",
        source="molinspect_policy",
        description="Named context() presets that expand to explicit radius, budget, and focus values.",
        parameters={
            "chemical_contacts": {
                "radius_A": 4.0,
                "budget": "small",
                "focus": ("contacts",),
            },
            "residue_environment": {
                "radius_A": 8.0,
                "budget": "medium",
                "focus": ("general",),
            },
            "ligand_binding_site": {
                "radius_A": 4.0,
                "budget": "medium",
                "focus": ("ligand_contact_shell",),
            },
            "metal_coordination": {
                "radius_A": 3.0,
                "budget": "small",
                "focus": ("metal_coordination",),
            },
            "protein_interface": {
                "radius_A": 5.0,
                "budget": "medium",
                "focus": ("interchain_interfaces",),
            },
            "broad_environment": {
                "radius_A": 12.0,
                "budget": "large",
                "focus": ("general",),
            },
        },
        limitations=("Context scales are retrieval presets, not structural-biology definitions.",),
    ),
    StructuralDefinition(
        id="distance_change_event",
        label="distance_change_event",
        kind="temporal_event",
        source="molinspect_policy",
        description="Reported local-context distance change between two compared frames.",
        parameters={"min_abs_delta_A": 0.5},
    ),
    StructuralDefinition(
        id="ligand_motion_event",
        label="ligand_motion_event",
        kind="temporal_event",
        source="molinspect_policy",
        description="Ligand, ion, contact-shell, or pocket distance change over adjacent sampled frames.",
        parameters={"min_abs_delta_A": 1.0},
    ),
    StructuralDefinition(
        id="selection_spread_event",
        label="selection_spread_event",
        kind="temporal_event",
        source="molinspect_policy",
        description="Selection radius-of-gyration spread increase or decrease over adjacent sampled frames.",
        parameters={"min_abs_delta_A": 0.5},
    ),
)

STRUCTURAL_DEFINITIONS: dict[str, StructuralDefinition] = {
    definition.id: definition for definition in _DEFINITIONS
}


def definition(definition_id: str) -> StructuralDefinition:
    """Return a structural definition by stable ID."""

    try:
        return STRUCTURAL_DEFINITIONS[definition_id]
    except KeyError as exc:
        raise KeyError(f"Unknown structural definition: {definition_id!r}") from exc


def definitions_by_kind(kind: DefinitionKind) -> tuple[StructuralDefinition, ...]:
    """Return definitions for one registry kind in stable order."""

    return tuple(definition for definition in _DEFINITIONS if definition.kind == kind)


def parameter(definition_id: str, parameter_name: str) -> Any:
    """Return a parameter from a structural definition."""

    definition_record = definition(definition_id)
    try:
        return definition_record.parameters[parameter_name]
    except KeyError as exc:
        raise KeyError(f"{definition_id!r} has no parameter {parameter_name!r}") from exc


def float_parameter(definition_id: str, parameter_name: str) -> float:
    """Return a definition parameter as a float."""

    return float(parameter(definition_id, parameter_name))


def int_parameter(definition_id: str, parameter_name: str) -> int:
    """Return a definition parameter as an int."""

    return int(parameter(definition_id, parameter_name))


def limitation(definition_id: str, index: int = 0) -> str:
    """Return one limitation string for a definition."""

    return definition(definition_id).limitations[index]


def relation_priority_for_type(relation_type: str) -> int:
    """Return stable sort priority for public relation labels."""

    definition_record = STRUCTURAL_DEFINITIONS.get(relation_type)
    if definition_record is None or definition_record.priority is None:
        return 99
    return definition_record.priority


def validate_definition_registry() -> None:
    """Validate internal definition/reference consistency."""

    if len(STRUCTURAL_DEFINITIONS) != len(_DEFINITIONS):
        raise ValueError("Structural definition IDs must be unique.")
    for definition_record in _DEFINITIONS:
        for reference_key in definition_record.reference_keys:
            if reference_key not in REFERENCE_CATALOG:
                raise ValueError(
                    f"{definition_record.id!r} references unknown key {reference_key!r}"
                )


validate_definition_registry()

CONTACT_CUTOFF_A = float_parameter("nonbonded_contact", "max_distance_A")
STERIC_CLASH_MIN_VDW_OVERLAP_A = float_parameter("steric_clash", "min_vdw_overlap_A")
INFERRED_BOND_VDW_FUDGE_FACTOR = float_parameter(
    "inferred_covalent_bond", "vdw_fudge_factor"
)
INFERRED_BOND_MIN_DISTANCE_A = float_parameter("inferred_covalent_bond", "minimum_distance_A")
PEPTIDE_BOND_C_N_MAX_DISTANCE_A = float_parameter(
    "inferred_covalent_bond", "peptide_c_n_max_distance_A"
)
METAL_COORDINATION_CUTOFF_A = float_parameter("metal_coordination", "max_distance_A")
HBOND_DISTANCE_CUTOFF_A = float_parameter("hydrogen_bond", "max_distance_A")
HBOND_MIN_ANGLE_DEG = float_parameter("hydrogen_bond", "min_angle_deg")
HBOND_DONOR_HYDROGEN_MAX_DISTANCE_A = float_parameter(
    "hydrogen_bond", "donor_hydrogen_max_distance_A"
)
TEMPORAL_HBOND_DONOR_ACCEPTOR_DISTANCE_A = float_parameter(
    "hydrogen_bond_occupancy", "donor_acceptor_distance_cutoff_A"
)
TEMPORAL_HBOND_DONOR_HYDROGEN_DISTANCE_A = float_parameter(
    "hydrogen_bond_occupancy", "donor_hydrogen_distance_cutoff_A"
)
TEMPORAL_HBOND_MIN_ANGLE_DEG = float_parameter(
    "hydrogen_bond_occupancy", "donor_hydrogen_acceptor_angle_min_deg"
)
LIGAND_STABILITY_STABLE_MIN_CONTACT_OCCUPANCY = float_parameter(
    "ligand_stability", "stable_min_contact_occupancy"
)
LIGAND_STABILITY_STABLE_MAX_RMSD_A = float_parameter(
    "ligand_stability", "stable_max_site_aligned_ligand_rmsd_A"
)
LIGAND_STABILITY_UNSTABLE_MAX_CONTACT_OCCUPANCY = float_parameter(
    "ligand_stability", "unstable_max_contact_occupancy"
)
LIGAND_STABILITY_UNSTABLE_MIN_RMSD_A = float_parameter(
    "ligand_stability", "unstable_min_site_aligned_ligand_rmsd_A"
)
STATE_MIN_RMSD_RANGE_A = float_parameter("conformational_states", "min_rmsd_range_A_for_two_states")
SALT_BRIDGE_CUTOFF_A = float_parameter("salt_bridge", "max_distance_A")
HYDROPHOBIC_CUTOFF_A = float_parameter("hydrophobic_contact", "max_distance_A")
PI_STACKING_CENTER_CUTOFF_A = float_parameter("pi_stacking", "max_center_distance_A")
PI_STACKING_MAX_ANGLE_DEVIATION_DEG = float_parameter(
    "pi_stacking", "max_angle_deviation_deg"
)
PI_STACKING_MAX_PARALLEL_OFFSET_A = float_parameter("pi_stacking", "max_parallel_offset_A")
WATER_BRIDGE_MIN_DISTANCE_A = float_parameter("water_bridge_candidate", "min_leg_distance_A")
WATER_BRIDGE_MAX_DISTANCE_A = float_parameter("water_bridge_candidate", "max_leg_distance_A")
EXPOSURE_NEIGHBOR_RADIUS_A = float_parameter("local_packing_density", "neighbor_radius_A")
EXPOSURE_DENSITY_LOWER_PERCENTILE = float_parameter(
    "local_packing_density", "lower_percentile"
)
EXPOSURE_DENSITY_UPPER_PERCENTILE = float_parameter(
    "local_packing_density", "upper_percentile"
)
FREESASA_BURIED_RELATIVE_SASA_MAX = float_parameter(
    "freesasa_exposure", "buried_relative_sasa_max"
)
FREESASA_EXPOSED_RELATIVE_SASA_MIN = float_parameter(
    "freesasa_exposure", "exposed_relative_sasa_min"
)
INTERFACE_DISTANCE_A = float_parameter("interchain_contact_interface", "contact_cutoff_A")
LIGAND_CONTACT_DISTANCE_A = float_parameter("ligand_contact_shell", "contact_cutoff_A")
MAX_RESOLVED_OBJECT_REFS = int_parameter("selection_region", "max_residue_objects")
DEFAULT_CONTEXT_RADIUS_A = float_parameter("context_radius", "default_radius_A")
BUDGET_LIMITS = {key: int(value) for key, value in definition("context_budget").parameters.items()}
DISTANCE_CHANGE_MIN_DELTA_A = float_parameter("distance_change_event", "min_abs_delta_A")
LIGAND_MOTION_MIN_DELTA_A = float_parameter("ligand_motion_event", "min_abs_delta_A")
SELECTION_SPREAD_MIN_DELTA_A = float_parameter("selection_spread_event", "min_abs_delta_A")
