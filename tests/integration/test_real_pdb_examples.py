from __future__ import annotations

import pytest

from molinspect import load
from tests.helpers import example_data_dir

DATA_DIR = example_data_dir()
REQUIRED = {
    "1crn.pdb": {"atoms": 327, "residues": 46, "chains": 1, "protein_residues": 46},
    "1ubq.pdb": {"atoms": 660, "residues": 134, "chains": 1, "protein_residues": 76},
    "4hhb.pdb": {"atoms": 4779, "residues": 801, "chains": 4, "protein_residues": 574},
}

pytestmark = pytest.mark.skipif(
    not all((DATA_DIR / filename).exists() for filename in REQUIRED),
    reason="real PDB example files are not present; run scripts/fetch_example_pdbs.py",
)


@pytest.mark.parametrize(("filename", "expected"), REQUIRED.items())
def test_real_pdb_examples_load_and_classify(filename, expected):
    session = load(structure=DATA_DIR / filename)
    summary = session.summary()

    assert summary.n_atoms == expected["atoms"]
    assert summary.n_residues == expected["residues"]
    assert summary.n_chains == expected["chains"]
    assert summary.n_frames == 1
    assert summary.mode == "static"
    assert session.objects(type="residue", limit=1000).count == expected["protein_residues"]


def test_known_residue_selection_on_real_pdbs():
    assert load(structure=DATA_DIR / "1crn.pdb").resolve_selection("A:1").resolved_objects == [
        "residue:A:1:THR"
    ]
    assert load(structure=DATA_DIR / "1ubq.pdb").resolve_selection("chain A and resid 1").resolved_objects == [
        "residue:A:1:MET"
    ]


def test_4hhb_chain_and_ligand_objects():
    session = load(structure=DATA_DIR / "4hhb.pdb")
    chains = session.objects(type="chain", limit=10)
    ligands = session.objects(type="ligand", limit=20)

    assert [obj.id for obj in chains.objects] == ["chain:A", "chain:B", "chain:C", "chain:D"]
    assert ligands.count == 6
    assert "ligand:A:142:HEM" in {obj.id for obj in ligands.objects}
    assert "ligand:B:147:PO4" in {obj.id for obj in ligands.objects}


def test_4hhb_heme_context_uses_real_structure():
    session = load(structure=DATA_DIR / "4hhb.pdb")

    located = session.locate("chain A and resname HEM")
    context = session.context("chain A and resname HEM", radius=4.0)

    assert located.resolved_objects == ["ligand:A:142:HEM"]
    assert located.location.near_ligands[0].object_id.startswith("ligand:")
    assert any(relation.type == "nonbonded_contact" for relation in context.relations)
    assert len(context.objects) < session.summary().n_residues


def test_4hhb_heme_iron_context_reports_metal_coordination():
    session = load(structure=DATA_DIR / "4hhb.pdb")

    located = session.locate("chain A and resid 87")
    context = session.context(
        "chain A and resname HEM and name FE",
        radius=3.0,
        focus="metal_coordination",
    )

    assert located.location.structural_profile
    assert located.location.local_packing
    if located.location.secondary_structure:
        assert located.location.secondary_structure
    assert "ligand:A:142:HEM" in located.location.ligand_contacts
    assert context.relations[0].type == "metal_coordination"
    assert context.relations[0].backend == "PLIP"
    assert context.relations[0].target == "residue:A:87:HIS"
    assert context.relations[0].source_atom == "A:142:HEM:FE"
    assert context.relations[0].target_atom == "A:87:HIS:NE2"


def test_4hhb_heme_iron_relation_timeline_uses_atom_pair_evidence():
    session = load(structure=DATA_DIR / "4hhb.pdb")

    result = session.timeline(
        metric="relation",
        selection1="chain A and resname HEM and name FE",
        selection2="chain A and resid 87",
    )

    assert result.summary["dominant_relation_type"] == "metal_coordination"
    assert result.sampled_values == [
        {
            "frame": 0,
            "relation_type": "metal_coordination",
            "distance_A": 2.143,
            "source_atom": "A:142:HEM:FE",
            "target_atom": "A:87:HIS:NE2",
        }
    ]


def test_real_pdb_relation_examples_cover_non_metal_relation_types():
    ubiquitin = load(structure=DATA_DIR / "1ubq.pdb")
    crambin = load(structure=DATA_DIR / "1crn.pdb")
    hemoglobin = load(structure=DATA_DIR / "4hhb.pdb")

    salt_bridge = ubiquitin.context("chain A and resid 11", radius=5.0, focus="salt_bridges").relations[0]
    hbond = crambin.context("chain A and resid 1", radius=5.0, focus="hydrogen_bonds").relations[0]
    hydrophobic = ubiquitin.context(
        "chain A and resid 3",
        radius=5.0,
        focus="hydrophobic_contacts",
    ).relations[0]
    aromatic = hemoglobin.context(
        "chain A and resid 33",
        radius=5.0,
        focus="pi_stacking",
    ).relations[0]

    assert (salt_bridge.type, salt_bridge.target, salt_bridge.source_atom, salt_bridge.target_atom) == (
        "salt_bridge",
        "residue:A:34:GLU",
        "A:11:LYS:NZ",
        "A:34:GLU:OE2",
    )
    assert hbond.type == "polar_contact_candidate"
    assert hbond.target == "residue:A:35:ILE"
    assert hbond.backend == "PDBe Arpeggio"
    assert hydrophobic.type == "hydrophobic_contact"
    assert hydrophobic.target == "residue:A:15:LEU"
    assert hydrophobic.backend == "PDBe Arpeggio"
    assert aromatic.type == "pi_stacking"
    assert aromatic.target == "residue:A:43:PHE"


def test_real_context_scales_prioritize_scientific_signal():
    ubiquitin = load(structure=DATA_DIR / "1ubq.pdb")
    hemoglobin = load(structure=DATA_DIR / "4hhb.pdb")

    lys_contacts = ubiquitin.context("chain A and resid 11", scale="chemical_contacts")
    his_environment = hemoglobin.context("chain A and resid 87", scale="residue_environment")
    chain_interface = hemoglobin.context("chain A", scale="protein_interface")

    assert lys_contacts.relations[0].type == "salt_bridge"
    assert lys_contacts.relations[0].target == "residue:A:34:GLU"
    assert "Top-ranked relation is salt_bridge" in lys_contacts.summary
    assert his_environment.relations[0].type == "metal_coordination"
    assert his_environment.relations[0].target == "ligand:A:142:HEM"
    assert "Top-ranked relation is metal_coordination" in his_environment.summary
    assert chain_interface.relations[0].target.startswith(
        ("residue:B:", "residue:C:", "residue:D:")
    )
    assert len(chain_interface.objects) < 80


def test_broad_real_selection_is_compacted():
    session = load(structure=DATA_DIR / "4hhb.pdb")

    resolved = session.resolve_selection("protein")
    context = session.context("chain A", radius=4.0)

    assert resolved.resolved_objects[0].startswith("selection_region:")
    assert "represented as one selection_region" in resolved.limitations[-1]
    assert context.objects[0].type == "selection_region"
    assert len(context.objects) < 50


def test_timeline_preserves_broad_selection_limitations():
    session = load(structure=DATA_DIR / "4hhb.pdb")

    result = session.timeline(metric="distance", selection1="protein", selection2="ligand")

    assert "represented as one selection_region" in " ".join(result.limitations)
