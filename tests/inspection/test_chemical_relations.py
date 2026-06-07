from __future__ import annotations

from types import SimpleNamespace
from typing import Any

import numpy as np

from molinspect import load
from molinspect.relations import relation_for_atomgroups


class FakeResidue:
    def __init__(self, resname: str, resid: int, chain: str = "A") -> None:
        self.resname = resname
        self.resid = resid
        self.segid = chain
        self.atoms: list[FakeAtom] = []


class FakeAtom:
    def __init__(
        self,
        name: str,
        element: str,
        residue: FakeResidue,
        position: tuple[float, float, float],
        chain: str = "A",
    ) -> None:
        self.name = name
        self.element = element
        self.residue = residue
        self.position = np.asarray(position, dtype=float)
        self.chainID = chain
        residue.atoms.append(self)


class FakeAtomGroup(list[FakeAtom]):
    @property
    def positions(self) -> np.ndarray:
        return np.asarray([atom.position for atom in self], dtype=float)


def _pair(atom_a: FakeAtom, atom_b: FakeAtom) -> Any:
    return SimpleNamespace(
        atom_a=atom_a,
        atom_b=atom_b,
        distance_A=round(float(np.linalg.norm(atom_a.position - atom_b.position)), 3),
    )


def _relation(
    atom_a: FakeAtom,
    atom_b: FakeAtom,
    source_atoms: FakeAtomGroup | None = None,
    target_atoms: FakeAtomGroup | None = None,
):
    source_group = source_atoms or FakeAtomGroup([atom_a])
    target_group = target_atoms or FakeAtomGroup([atom_b])
    return relation_for_atomgroups(source_group, target_group, "source", "target", _pair(atom_a, atom_b))


def test_metal_coordination_relation_has_specific_cutoff_and_limitation():
    heme = FakeResidue("HEM", 142)
    histidine = FakeResidue("HIS", 87)
    iron = FakeAtom("FE", "FE", heme, (0.0, 0.0, 0.0))
    nitrogen = FakeAtom("NE2", "N", histidine, (2.1, 0.0, 0.0))

    relation = _relation(iron, nitrogen)

    assert relation.type == "metal_coordination"
    assert relation.category == "coordination"
    assert relation.confidence == "geometry"
    assert relation.cutoff_A == 3.0
    assert relation.min_distance_A == 2.1
    assert "not full valence geometry" in relation.limitations[0]


def test_salt_bridge_uses_charged_atom_names_not_only_residue_class():
    lysine = FakeResidue("LYS", 10)
    aspartate = FakeResidue("ASP", 30)
    nz = FakeAtom("NZ", "N", lysine, (0.0, 0.0, 0.0))
    od1 = FakeAtom("OD1", "O", aspartate, (3.2, 0.0, 0.0))

    relation = _relation(nz, od1)

    assert relation.type == "salt_bridge"
    assert relation.category == "electrostatic"
    assert relation.cutoff_A == 5.5


def test_hydrogen_bond_requires_explicit_angle_when_hydrogen_is_available():
    donor_residue = FakeResidue("ASN", 20)
    acceptor_residue = FakeResidue("ASP", 45)
    donor = FakeAtom("ND2", "N", donor_residue, (0.0, 0.0, 0.0))
    hydrogen = FakeAtom("HD21", "H", donor_residue, (1.0, 0.0, 0.0))
    acceptor = FakeAtom("OD1", "O", acceptor_residue, (2.8, 0.0, 0.0))

    relation = _relation(
        donor,
        acceptor,
        source_atoms=FakeAtomGroup([donor, hydrogen]),
        target_atoms=FakeAtomGroup([acceptor]),
    )

    assert relation.type == "hydrogen_bond"
    assert relation.category == "polar"
    assert relation.cutoff_A == 4.1
    assert relation.angle_deg == 180.0
    assert relation.evidence[1].metric == "donor_hydrogen_acceptor_angle"


def test_polar_contact_candidate_is_reserved_for_missing_angle_validation():
    donor_residue = FakeResidue("ASN", 20)
    acceptor_residue = FakeResidue("ASP", 45)
    donor = FakeAtom("ND2", "N", donor_residue, (0.0, 0.0, 0.0))
    acceptor = FakeAtom("OD1", "O", acceptor_residue, (2.8, 0.0, 0.0))

    relation = _relation(donor, acceptor)

    assert relation.type == "polar_contact_candidate"
    assert relation.confidence == "candidate"
    assert relation.angle_deg is None
    assert "lack explicit donor/acceptor and angle validation" in relation.limitations[0]


def test_aromatic_and_hydrophobic_contacts_are_separate_from_generic_contacts():
    phenylalanine = FakeResidue("PHE", 10)
    tyrosine = FakeResidue("TYR", 40)
    leu = FakeResidue("LEU", 50)
    val = FakeResidue("VAL", 75)
    phe_ring = _add_planar_ring(phenylalanine, z=0.0)[0]
    tyr_ring = _add_planar_ring(tyrosine, z=3.5)[0]
    leu_carbon = FakeAtom("CD1", "C", leu, (0.0, 0.0, 0.0))
    val_carbon = FakeAtom("CG1", "C", val, (3.8, 0.0, 0.0))

    aromatic = _relation(phe_ring, tyr_ring)
    hydrophobic = _relation(leu_carbon, val_carbon)

    assert aromatic.type == "pi_stacking"
    assert aromatic.cutoff_A == 5.5
    assert hydrophobic.type == "hydrophobic_contact"
    assert hydrophobic.cutoff_A == 4.0


def test_steric_clash_excludes_adjacent_backbone_peptide_bonds():
    residue_a = FakeResidue("ALA", 1)
    residue_b = FakeResidue("GLY", 2)
    residue_c = FakeResidue("SER", 9)
    backbone_c = FakeAtom("C", "C", residue_a, (0.0, 0.0, 0.0))
    backbone_n = FakeAtom("N", "N", residue_b, (1.5, 0.0, 0.0))
    nonbonded_c = FakeAtom("CB", "C", residue_c, (0.0, 1.6, 0.0))

    peptide_bond = _relation(backbone_c, backbone_n)
    clash = _relation(backbone_c, nonbonded_c)

    assert peptide_bond.type == "inferred_covalent_bond"
    assert peptide_bond.backend == "molinspect_heuristic"
    assert clash.type == "steric_clash"
    assert clash.category == "steric"
    assert clash.cutoff_A == 0.4


def test_near_relation_is_explicit_when_pair_is_outside_contact_cutoff():
    residue_a = FakeResidue("ALA", 1)
    residue_b = FakeResidue("GLY", 10)
    atom_a = FakeAtom("CB", "C", residue_a, (0.0, 0.0, 0.0))
    atom_b = FakeAtom("CA", "C", residue_b, (6.0, 0.0, 0.0))

    relation = _relation(atom_a, atom_b)

    assert relation.type == "near"
    assert relation.category == "proximity"
    assert relation.min_distance_A == 6.0


def test_context_reports_water_bridge_candidate(tmp_path):
    pdb = tmp_path / "water_bridge.pdb"
    pdb.write_text(
        _pdb_line("ATOM", 1, "N", "ASN", "A", 1, 0.0, 0.0, 0.0, "N")
        + _pdb_line("ATOM", 2, "O", "ASP", "A", 2, 5.6, 0.0, 0.0, "O")
        + _pdb_line("HETATM", 3, "O", "HOH", "A", 100, 2.8, 0.0, 0.0, "O")
        + "END\n"
    )
    session = load(structure=pdb)

    result = session.context(
        "chain A and resid 1",
        radius=6.0,
        focus="water_bridges",
    )

    relation = result.relations[0]
    assert relation.type == "water_bridge_candidate"
    assert relation.target == "residue:A:2:ASP"
    assert relation.category == "water_mediated"
    assert relation.confidence == "candidate"
    assert relation.cutoff_A == 4.1
    assert relation.evidence[1].metric == "source_to_water_distance"
    assert relation.evidence[2].metric == "water_to_target_distance"
    assert relation.evidence[3].value == "water:A:100:HOH"


def _add_planar_ring(residue: FakeResidue, z: float) -> list[FakeAtom]:
    names = ("CG", "CD1", "CE1", "CZ", "CE2", "CD2")
    coordinates = (
        (1.4, 0.0, z),
        (0.7, 1.212, z),
        (-0.7, 1.212, z),
        (-1.4, 0.0, z),
        (-0.7, -1.212, z),
        (0.7, -1.212, z),
    )
    return [
        FakeAtom(name, "C", residue, coordinates[index])
        for index, name in enumerate(names)
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
