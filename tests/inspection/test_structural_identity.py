from __future__ import annotations

import sys
from pathlib import Path
from types import SimpleNamespace

from molinspect import load
from molinspect.objects import _object_ref_for_atom, object_id_for_atom, object_id_for_residue
from tests.helpers import tiny_pdb_text


class FakeAtoms(list):
    @property
    def icodes(self):
        return [getattr(atom, "icode", "") for atom in self]


class FakeResidue:
    def __init__(
        self,
        *,
        resname: str = "GLY",
        resid: int = 57,
        icode: str = "",
        chain_id: str = "",
        segid: str = "",
    ) -> None:
        self.resname = resname
        self.resid = resid
        self.icode = icode
        self.segid = segid
        self.atoms = FakeAtoms()
        atom = FakeAtom(
            residue=self,
            name="CA",
            ix=4,
            chain_id=chain_id,
            segid=segid,
            altloc="B",
            icode=icode,
        )
        self.atoms.append(atom)


class FakeAtom:
    def __init__(
        self,
        *,
        residue: FakeResidue,
        name: str,
        ix: int,
        chain_id: str = "",
        segid: str = "",
        altloc: str = "",
        icode: str = "",
    ) -> None:
        self.residue = residue
        self.name = name
        self.ix = ix
        self.chainID = chain_id
        self.segid = segid
        self.altLoc = altloc
        self.icode = icode


def test_object_ids_include_insertion_code_and_fallback_segment():
    residue = FakeResidue(chain_id="", segid="LONGSEG", icode="A")
    atom = residue.atoms[0]

    assert object_id_for_residue(residue, object_type="residue") == "residue:LONGSEG:57A:GLY"
    assert object_id_for_atom(atom) == "atom:LONGSEG:57A:GLY:CA:4"


def test_atom_object_ref_exposes_altloc_and_icode():
    residue = FakeResidue(chain_id="AB", segid="", icode="C")
    atom = residue.atoms[0]

    ref = _object_ref_for_atom(atom)

    assert ref.id == "atom:AB:57C:GLY:CA:4"
    assert ref.chain == "AB"
    assert ref.icode == "C"
    assert ref.altloc == "B"


def test_freesasa_profile_is_used_when_available(monkeypatch, tiny_pdb):
    class FakeArea:
        total = 65.0
        relativeTotal = 42.0

    class FakeResult:
        def residueAreas(self):
            return {"A": {"1": FakeArea()}}

    fake_freesasa = SimpleNamespace(
        Structure=lambda path: object(),
        calc=lambda structure: FakeResult(),
    )
    monkeypatch.setitem(sys.modules, "freesasa", fake_freesasa)

    result = load(structure=tiny_pdb).locate("chain A and resid 1")
    profile = result.location.structural_profile

    assert profile is not None
    assert profile.exposure_source == "freesasa"
    assert profile.sasa_A2 == 65.0
    assert profile.relative_sasa == 0.42
    assert result.location.exposure_status == "exposed"
    assert result.location.surface_status == "surface"
    assert "FreeSASA exact SASA is unavailable" not in " ".join(result.limitations)


def test_mmcif_loader_uses_optional_gemmi_adapter(monkeypatch, tmp_path):
    mmcif_path = tmp_path / "tiny.cif"
    mmcif_path.write_text("data_tiny\n#\n")

    class FakeGemmiStructure:
        def write_pdb(self, output_path: str) -> None:
            Path(output_path).write_text(tiny_pdb_text())

    fake_gemmi = SimpleNamespace(read_structure=lambda path: FakeGemmiStructure())
    monkeypatch.setitem(sys.modules, "gemmi", fake_gemmi)

    session = load(structure=mmcif_path)

    assert session.summary().n_atoms == 11
    assert session.resolve_selection("A:1").resolved_objects == ["residue:A:1:ALA"]
    assert session.objects(type="ligand").objects[0].id == "ligand:A:101:ATP"
