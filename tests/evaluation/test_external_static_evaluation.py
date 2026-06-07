from __future__ import annotations

from molinspect import load
from tests.helpers_static_spatial import (
    ATOM3D_PIP_INTERFACE_CUTOFF_A,
    atom3d_pip_interface_residues,
    context_residue_keys,
    external_static_benchmark_sources,
    parse_biolip_record,
    parse_p2rank_dataset_items,
    residue_key,
    score_pocket_center_hit,
    score_residue_retrieval,
)
from molinspect.schemas import ContextResult, ObjectRef, Relation


def test_external_static_benchmark_registry_is_static_and_not_private():
    sources = external_static_benchmark_sources()
    ids = {source.id for source in sources}

    assert {"biolip2", "cameo_ligand_binding_site", "p2rank_datasets", "atom3d_pip"}.issubset(ids)
    assert all("temporal" not in source.fit_for_molinspect.lower() for source in sources)
    assert all(source.url.startswith("https://") for source in sources)
    assert any(source.priority == "primary" for source in sources)


def test_residue_key_normalizes_object_ids_and_chain_resid_labels():
    assert residue_key("residue:A:87:HIS").label == "A:87"
    assert residue_key("A:87").label == "A:87"
    assert residue_key("B:12:GLY").resname == "GLY"


def test_context_residue_keys_extracts_objects_and_relation_endpoints():
    context = ContextResult(
        selection="resname HEM",
        frame=0,
        radius_A=4.0,
        focus=["ligand_contact_shell"],
        objects=[
            ObjectRef(id="ligand:A:142:HEM", type="ligand", name="HEM"),
            ObjectRef(id="residue:A:87:HIS", type="residue", name="HIS"),
        ],
        relations=[
            Relation(
                source="ligand:A:142:HEM",
                target="residue:A:92:LEU",
                type="nonbonded_contact",
                min_distance_A=3.5,
            )
        ],
        summary="2 residue keys retrieved.",
    )

    assert {key.label for key in context_residue_keys(context)} == {"A:87", "A:92"}


def test_score_residue_retrieval_reports_exact_error_sets():
    score = score_residue_retrieval(
        predicted=["residue:A:87:HIS", "A:92", "A:93"],
        gold=["A:87", "residue:A:92:LEU", "A:94"],
    )

    assert score.true_positives == ["A:87", "A:92"]
    assert score.false_positives == ["A:93"]
    assert score.false_negatives == ["A:94"]
    assert score.precision == 0.666667
    assert score.recall == 0.666667
    assert score.f1 == 0.666667


def test_parse_biolip_record_uses_pdb_numbered_binding_site_column():
    line = (
        "966c\tA\t1.90\tBS06\tRS2\tA\t1\t"
        "N180 L181 A182\tN73 L74 A75\tM236\tE219\t3.4.24.-\t0004222\t"
        "ki=23nM (RS2)\tKi=23nM (RS2)\t\t\tP03956\t10074939\t236\tSEQUENCE"
    )

    record = parse_biolip_record(line)

    assert record.pdb_id == "966c"
    assert record.receptor_chain == "A"
    assert record.resolution_A == 1.9
    assert record.binding_site_code == "BS06"
    assert record.ligand_auth_resid == "236"
    assert record.ligand_selection == "chain A and resid 236 and resname RS2"
    assert [key.label for key in record.binding_site_residues] == ["A:180", "A:181", "A:182"]


def test_parse_p2rank_dataset_items_accepts_comments_and_whitespace():
    text = """
    # subset
    chen11/1abcA.pdb chen11/2defB.pdb

    holo4k/3ghi.pdb  # inline comment
    """

    assert parse_p2rank_dataset_items(text) == [
        "chen11/1abcA.pdb",
        "chen11/2defB.pdb",
        "holo4k/3ghi.pdb",
    ]


def test_score_pocket_center_hit_uses_angstrom_distance_cutoff():
    hit = score_pocket_center_hit((0, 0, 0), (1, 2, 2), cutoff_A=3.0)
    miss = score_pocket_center_hit((0, 0, 0), (0, 0, 4.1), cutoff_A=4.0)

    assert hit.distance_A == 3.0
    assert hit.hit is True
    assert miss.distance_A == 4.1
    assert miss.hit is False


def test_atom3d_pip_interface_residues_use_six_angstrom_heavy_atom_rule(tmp_path):
    pdb = tmp_path / "interface.pdb"
    pdb.write_text(
        _pdb_line("ATOM", 1, "CA", "ALA", "A", 1, 0.0, 0.0, 0.0, "C")
        + _pdb_line("ATOM", 2, "CB", "ALA", "A", 1, 1.0, 0.0, 0.0, "C")
        + _pdb_line("ATOM", 3, "CA", "GLY", "B", 5, 6.5, 0.0, 0.0, "C")
        + _pdb_line("ATOM", 4, "CA", "GLY", "B", 6, 10.5, 0.0, 0.0, "C")
        + "END\n"
    )
    session = load(structure=pdb)

    assert ATOM3D_PIP_INTERFACE_CUTOFF_A == 6.0
    assert [key.label for key in atom3d_pip_interface_residues(session.universe, "A", "B")] == [
        "B:5"
    ]
    assert [key.label for key in atom3d_pip_interface_residues(session.universe, "B", "A")] == [
        "A:1"
    ]


def test_atom3d_pip_interface_residues_preserve_insertion_codes(tmp_path):
    pdb = tmp_path / "insertion_code_interface.pdb"
    pdb.write_text(
        _pdb_line("ATOM", 1, "CA", "ALA", "A", 221, 0.0, 0.0, 0.0, "C", icode="A")
        + _pdb_line("ATOM", 2, "CB", "ALA", "A", 221, 1.0, 0.0, 0.0, "C", icode="A")
        + _pdb_line("ATOM", 3, "CA", "GLY", "B", 5, 6.5, 0.0, 0.0, "C")
        + "END\n"
    )
    session = load(structure=pdb)

    predicted = [
        key for key in context_residue_keys(session.context("protein and chain B", radius=6.0))
        if key.chain == "A"
    ]
    gold = atom3d_pip_interface_residues(session.universe, "B", "A")
    score = score_residue_retrieval(predicted, gold)

    assert [key.label for key in gold] == ["A:221A"]
    assert "A:221A" in {key.label for key in predicted}
    assert score.false_positives == []
    assert score.false_negatives == []


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
    icode: str = " ",
) -> str:
    return (
        f"{record:<6}{serial:5d} {name:^4} {resname:>3} {chain:1}{resid:4d}{icode:1}   "
        f"{x:8.3f}{y:8.3f}{z:8.3f}{1.00:6.2f}{20.00:6.2f}          {element:>2}\n"
    )
