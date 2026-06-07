from __future__ import annotations

import pytest

from molinspect import load
from molinspect.errors import SelectionResolutionError


def test_chain_selection_is_translated(tiny_pdb):
    session = load(structure=tiny_pdb)

    resolved = session.resolve_selection("chain A and resid 1")

    assert resolved.expression == "chainID A and resid 1"
    assert resolved.n_atoms == 3
    assert resolved.resolved_objects == ["residue:A:1:ALA"]


def test_short_residue_alias(tiny_pdb):
    session = load(structure=tiny_pdb)

    resolved = session.resolve_selection("A:1")

    assert resolved.expression == "chainID A and resid 1"
    assert resolved.resolved_objects == ["residue:A:1:ALA"]


def test_residue_range_selection_is_public_notation(tiny_pdb):
    session = load(structure=tiny_pdb)

    resolved = session.resolve_selection("chain A and resid 1-2")

    assert resolved.expression == "chainID A and resid 1-2"
    assert resolved.resolved_objects == ["residue:A:1:ALA", "residue:A:2:GLY"]


def test_ligand_water_and_ion_aliases(tiny_pdb):
    session = load(structure=tiny_pdb)

    ligand = session.resolve_selection("ligand")
    water = session.resolve_selection("water")
    ion = session.resolve_selection("ion")

    assert ligand.n_atoms == 3
    assert ligand.resolved_objects == ["ligand:A:101:ATP"]
    assert water.n_atoms == 1
    assert water.resolved_objects == ["water:A:201:HOH"]
    assert ion.n_atoms == 1
    assert ion.resolved_objects == ["ion:A:301:NA"]


def test_returned_object_ids_can_be_selected_again(tiny_pdb):
    session = load(structure=tiny_pdb)

    residue = session.resolve_selection("residue:A:1:ALA")
    ligand = session.resolve_selection("ligand:A:101:ATP")
    ion = session.resolve_selection("ion:A:301:NA")
    water = session.resolve_selection("water:A:201:HOH")
    atom = session.resolve_selection("atom:A:101:ATP:P:6")
    chain = session.resolve_selection("chain:A")

    assert residue.expression == "chainID A and resid 1 and resname ALA"
    assert residue.resolved_objects == ["residue:A:1:ALA"]
    assert ligand.resolved_objects == ["ligand:A:101:ATP"]
    assert ion.resolved_objects == ["ion:A:301:NA"]
    assert water.resolved_objects == ["water:A:201:HOH"]
    assert atom.resolved_objects == ["atom:A:101:ATP:P:6"]
    assert chain.resolved_objects == [
        "residue:A:1:ALA",
        "residue:A:2:GLY",
        "ligand:A:101:ATP",
        "water:A:201:HOH",
        "ion:A:301:NA",
    ]


def test_around_of_alias_is_translated(tiny_pdb):
    session = load(structure=tiny_pdb)

    resolved = session.resolve_selection("around 3 of chain A and resid 2")

    assert resolved.expression == "around 3 (chainID A and resid 2)"
    assert resolved.n_atoms == 3
    assert resolved.resolved_objects == ["residue:A:1:ALA", "ligand:A:101:ATP"]


def test_empty_selection_raises_clear_error(tiny_pdb):
    session = load(structure=tiny_pdb)

    with pytest.raises(SelectionResolutionError, match="resolved to zero atoms"):
        session.resolve_selection("chain Z")


def test_selection_error_explains_public_notation(tiny_pdb):
    session = load(structure=tiny_pdb)

    with pytest.raises(SelectionResolutionError, match="chain A and resid 57"):
        session.resolve_selection("")
