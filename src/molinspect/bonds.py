"""Covalent-bond provenance and steric-overlap helpers."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Literal

from .definitions import (
    INFERRED_BOND_MIN_DISTANCE_A,
    INFERRED_BOND_VDW_FUDGE_FACTOR,
    PEPTIDE_BOND_C_N_MAX_DISTANCE_A,
)

BondSource = Literal["topology_bond", "inferred_covalent_bond"]


@dataclass(frozen=True, slots=True)
class BondProvenance:
    """Covalent-bond evidence for one atom pair."""

    definition_id: BondSource
    confidence: Literal["topology", "inferred"]
    method: str
    cutoff_A: float | None = None


VDW_RADII_A: dict[str, float] = {
    "H": 1.20,
    "C": 1.70,
    "N": 1.55,
    "O": 1.52,
    "F": 1.47,
    "P": 1.80,
    "S": 1.80,
    "CL": 1.75,
    "BR": 1.85,
    "I": 1.98,
    "FE": 2.00,
    "ZN": 2.10,
    "MG": 1.73,
    "MN": 2.05,
    "CA": 2.31,
    "CU": 1.40,
    "CO": 2.00,
    "NI": 1.63,
    "NA": 2.27,
    "K": 2.75,
}


def atom_element(atom: Any) -> str:
    """Return a best-effort element symbol for atom-like objects."""

    element = str(getattr(atom, "element", "")).strip().upper()
    if element:
        return element
    name = str(getattr(atom, "name", "")).strip().upper()
    letters = "".join(character for character in name if character.isalpha())
    if len(letters) >= 2 and letters[:2] in VDW_RADII_A:
        return letters[:2]
    return letters[:1]


def bond_provenance_for_pair(pair: Any) -> BondProvenance | None:
    """Return bond provenance when the atom pair is covalently connected."""

    if _topology_has_bond(pair.atom_a, pair.atom_b):
        return BondProvenance(definition_id="topology_bond", confidence="topology", method="topology")

    peptide_cutoff = _peptide_bond_cutoff_A(pair)
    if peptide_cutoff is not None:
        return BondProvenance(
            definition_id="inferred_covalent_bond",
            confidence="inferred",
            method="polymer_backbone_c_n",
            cutoff_A=peptide_cutoff,
        )

    inferred_cutoff = _same_residue_inferred_bond_cutoff_A(pair)
    if inferred_cutoff is not None:
        return BondProvenance(
            definition_id="inferred_covalent_bond",
            confidence="inferred",
            method="same_residue_vdw_fudge",
            cutoff_A=inferred_cutoff,
        )
    return None


def vdw_overlap_A(atom_a: Any, atom_b: Any, distance_A: float) -> float | None:
    """Return non-bonded van der Waals overlap in angstroms."""

    radius_a = VDW_RADII_A.get(atom_element(atom_a))
    radius_b = VDW_RADII_A.get(atom_element(atom_b))
    if radius_a is None or radius_b is None:
        return None
    return round(radius_a + radius_b - distance_A, 3)


def _topology_has_bond(atom_a: Any, atom_b: Any) -> bool:
    try:
        return atom_b in atom_a.bonded_atoms
    except Exception:
        pass

    ix_a = getattr(atom_a, "ix", None)
    ix_b = getattr(atom_b, "ix", None)
    universe = getattr(atom_a, "universe", None) or getattr(atom_b, "universe", None)
    if ix_a is None or ix_b is None or universe is None:
        return False

    try:
        bond_indices = universe.bonds.indices
    except Exception:
        return False
    pair = {int(ix_a), int(ix_b)}
    return any({int(first), int(second)} == pair for first, second in bond_indices)


def _peptide_bond_cutoff_A(pair: Any) -> float | None:
    residue_a = pair.atom_a.residue
    residue_b = pair.atom_b.residue
    if residue_a is residue_b:
        return None
    names = {str(pair.atom_a.name).strip().upper(), str(pair.atom_b.name).strip().upper()}
    if names != {"C", "N"}:
        return None
    if _chain_for_residue(residue_a) != _chain_for_residue(residue_b):
        return None
    try:
        adjacent = abs(int(residue_a.resid) - int(residue_b.resid)) == 1
    except (TypeError, ValueError):
        return None
    if adjacent and pair.distance_A <= PEPTIDE_BOND_C_N_MAX_DISTANCE_A:
        return PEPTIDE_BOND_C_N_MAX_DISTANCE_A
    return None


def _same_residue_inferred_bond_cutoff_A(pair: Any) -> float | None:
    if pair.atom_a.residue is not pair.atom_b.residue:
        return None
    if pair.distance_A < INFERRED_BOND_MIN_DISTANCE_A:
        return None
    radius_a = VDW_RADII_A.get(atom_element(pair.atom_a))
    radius_b = VDW_RADII_A.get(atom_element(pair.atom_b))
    if radius_a is None or radius_b is None:
        return None
    cutoff = round(INFERRED_BOND_VDW_FUDGE_FACTOR * (radius_a + radius_b), 3)
    return cutoff if pair.distance_A <= cutoff else None


def _chain_for_residue(residue: Any) -> str:
    for atom in residue.atoms:
        chain = str(getattr(atom, "chainID", "")).strip()
        if chain:
            return chain
    segment = str(getattr(residue, "segid", "")).strip()
    return segment or "_"
