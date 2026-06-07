from __future__ import annotations

import pytest

from molinspect import load
from molinspect.annotations import AnnotationStore, ResidueAnnotation
from molinspect.objects import object_id_for_residue
from tests.helpers import example_data_dir

DATA_DIR = example_data_dir()
PDB_4HHB = DATA_DIR / "4hhb.pdb"


pytestmark_real = pytest.mark.skipif(
    not PDB_4HHB.exists(),
    reason="4HHB example file is not present; run scripts/fetch_example_pdbs.py",
)


def test_secondary_structure_strand_and_loop_landmarks_are_discoverable_and_selectable(tmp_path):
    pdb = tmp_path / "three_residues.pdb"
    pdb.write_text(
        _pdb_line("ATOM", 1, "CA", "ALA", "A", 1, 0.0, 0.0, 0.0, "C")
        + _pdb_line("ATOM", 2, "CA", "VAL", "A", 2, 3.8, 0.0, 0.0, "C")
        + _pdb_line("ATOM", 3, "CA", "GLY", "A", 3, 7.6, 0.0, 0.0, "C")
        + "END\n"
    )
    session = load(structure=pdb)
    residues = list(session.universe.residues)
    session.world._annotations = AnnotationStore(
        by_residue_ix={
            residues[0].ix: ResidueAnnotation(
                object_id=object_id_for_residue(residues[0]),
                secondary_structure="alpha_helix",
                secondary_structure_code="H",
                secondary_structure_source="test_dssp",
            ),
            residues[1].ix: ResidueAnnotation(
                object_id=object_id_for_residue(residues[1]),
                secondary_structure="beta_strand",
                secondary_structure_code="E",
                secondary_structure_source="test_dssp",
            ),
            residues[2].ix: ResidueAnnotation(
                object_id=object_id_for_residue(residues[2]),
                secondary_structure="loop",
                secondary_structure_code=" ",
                secondary_structure_source="test_dssp",
            ),
        },
        limitations=(),
    )
    session.world._landmarks = None

    landmarks = session.objects(type=["secondary_structure", "loop"], limit=10)
    landmark_ids = {obj.id for obj in landmarks.objects}
    helix_id = "secondary_structure:A:alpha_helix:1-1"
    strand_id = "secondary_structure:A:beta_strand:2-2"
    loop_id = "loop:A:3-3"

    assert landmark_ids == {helix_id, strand_id, loop_id}
    assert session.resolve_selection(helix_id).resolved_objects == ["residue:A:1:ALA"]
    assert session.resolve_selection(strand_id).resolved_objects == ["residue:A:2:VAL"]
    located = session.locate("chain A and resid 3", include_metrics=False)
    assert located.location.secondary_structure_element == loop_id
    assert loop_id in located.location.landmark_memberships


def test_broad_ligand_contact_shell_selection_compacts_to_selection_region(tmp_path):
    pdb = tmp_path / "broad_site.pdb"
    lines = []
    serial = 1
    for resid in range(1, 31):
        lines.append(_pdb_line("ATOM", serial, "CA", "ALA", "A", resid, float(resid), 0.0, 0.0, "C"))
        serial += 1
    lines.append(_pdb_line("HETATM", serial, "P", "ATP", "A", 900, 40.0, 0.0, 0.0, "P"))
    pdb.write_text("".join(lines) + "END\n")
    session = load(structure=pdb)
    residues = list(session.universe.residues)
    ligand_id = "ligand:A:900:ATP"
    session.world._annotations = AnnotationStore(
        by_residue_ix={
            residue.ix: ResidueAnnotation(
                object_id=object_id_for_residue(residue),
                ligand_contact_ids=(ligand_id,) if int(residue.resid) <= 30 else (),
            )
            for residue in residues
        },
        limitations=(),
    )
    session.world._landmarks = None

    shell = session.objects(type="ligand_contact_shell", contains="ATP").objects[0]
    resolved = session.resolve_selection(shell.id)

    assert shell.annotations["lining_residue_count"] == 30
    assert resolved.resolved_objects[0].startswith("selection_region:")
    assert "represented as one selection_region" in " ".join(resolved.limitations)


@pytestmark_real
def test_4hhb_ligand_contact_shells_and_interchain_interfaces_are_first_class_objects():
    session = load(structure=PDB_4HHB)

    ligand_shells = session.objects(type="ligand_contact_shell", contains="HEM", limit=10)
    pockets = session.objects(type="pocket", contains="HEM", limit=10)
    interfaces = session.objects(type="interchain_contact_interface", limit=10)

    assert "ligand_contact_shell:A:142:HEM" in {obj.id for obj in ligand_shells.objects}
    assert pockets.count == 0
    assert interfaces.count >= 1
    assert all(obj.annotations["method"] for obj in ligand_shells.objects + interfaces.objects)


@pytestmark_real
def test_4hhb_ligand_contact_shell_membership_is_reflected_in_locate_and_context():
    session = load(structure=PDB_4HHB)

    located = session.locate("chain A and resid 87")
    context = session.context("chain A and resid 87", radius=4.0, focus="ligand_contact_shell")
    heme_shell = session.resolve_selection("ligand_contact_shell:A:142:HEM")

    assert located.location.ligand_contact_shell_ids == ["ligand_contact_shell:A:142:HEM"]
    assert located.location.pocket_ids == []
    assert "ligand_contact_shell:A:142:HEM" in located.location.landmark_memberships
    assert context.relations[0].target == "ligand:A:142:HEM"
    assert context.objects[0].annotations["ligand_contact_shell_ids"] == [
        "ligand_contact_shell:A:142:HEM"
    ]
    assert "ligand:A:142:HEM" in heme_shell.resolved_objects


@pytestmark_real
def test_4hhb_interface_focus_returns_interchain_landmark_aware_context():
    session = load(structure=PDB_4HHB)
    located = session.locate("chain A and resid 92")
    context = session.context("chain A and resid 92", radius=5.0, focus="interchain_interfaces")

    assert located.location.interchain_contact_interface_ids == [
        "interchain_contact_interface:A-D"
    ]
    assert context.relations[0].target.startswith("residue:D:")
    assert context.objects[0].annotations["interchain_contact_interface_ids"] == [
        "interchain_contact_interface:A-D"
    ]


def _pdb_line(
    record: str,
    serial: int,
    name: str,
    resname: str,
    chain: str,
    resid: int,
    x: float,
    y: float,
    z: float,
    element: str,
) -> str:
    return (
        f"{record:<6}{serial:5d} {name:^4} {resname:>3} {chain:1}{resid:4d}    "
        f"{x:8.3f}{y:8.3f}{z:8.3f}{1.00:6.2f}{20.00:6.2f}          {element:>2}\n"
    )
