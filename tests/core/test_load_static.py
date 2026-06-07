from __future__ import annotations

from molinspect import InspectionSession, load


def test_load_static_pdb(tiny_pdb):
    session = load(structure=tiny_pdb)

    assert isinstance(session, InspectionSession)
    assert session.n_frames == 1
    assert session.summary().n_atoms == 11
    assert session.summary().n_residues == 5
    assert session.summary().n_chains == 1
    assert session.summary().mode == "static"
    annotations = session.summary().available_annotations
    assert "basic_topology" in annotations
    assert "local_packing_density" in annotations
    assert "ligand_contact_shell" in annotations
    assert "freesasa_exposure" in annotations
    assert session.target_id.startswith("sha256:")
