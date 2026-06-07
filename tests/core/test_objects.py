from __future__ import annotations

import pytest

from molinspect import ObjectQueryError, load


def test_default_objects_omit_atoms(tiny_pdb):
    session = load(structure=tiny_pdb)

    result = session.objects()

    assert result.count == 6
    assert result.returned == 6
    assert not result.truncated
    assert "Atom objects are omitted by default" in result.limitations[0]
    assert {obj.type for obj in result.objects} == {"chain", "residue", "ligand", "water", "ion"}


def test_residue_objects(tiny_pdb):
    session = load(structure=tiny_pdb)

    result = session.objects(type="residue")

    assert result.count == 2
    assert [obj.id for obj in result.objects] == ["residue:A:1:ALA", "residue:A:2:GLY"]


def test_ligand_and_ion_objects(tiny_pdb):
    session = load(structure=tiny_pdb)

    result = session.objects(type=["ligand", "ion"])

    assert [obj.id for obj in result.objects] == ["ligand:A:101:ATP", "ion:A:301:NA"]
    assert [obj.atom_count for obj in result.objects] == [3, 1]


def test_atom_objects_are_explicit_and_limited(tiny_pdb):
    session = load(structure=tiny_pdb)

    result = session.objects(type="atom", limit=2)

    assert result.count == 11
    assert result.returned == 2
    assert result.truncated
    assert result.objects[0].id == "atom:A:1:ALA:N:0"


def test_object_contains_literal_filter(tiny_pdb):
    session = load(structure=tiny_pdb)

    result = session.objects(contains="ATP")

    assert result.count == 1
    assert result.objects[0].id == "ligand:A:101:ATP"


def test_object_contains_rejects_empty_literal(tiny_pdb):
    session = load(structure=tiny_pdb)

    with pytest.raises(ObjectQueryError, match="non-empty literal substring"):
        session.objects(contains=" ")


def test_object_type_error_explains_supported_notation(tiny_pdb):
    session = load(structure=tiny_pdb)

    with pytest.raises(ObjectQueryError, match="Supported object types"):
        session.objects(type="domain")
