"""Shared chemical relation semantics for context, timeline, and compare."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any

import numpy as np

from .bonds import atom_element, bond_provenance_for_pair, vdw_overlap_A
from .schemas import EvidenceItem, Relation
from .definitions import (
    AROMATIC_RING_ATOMS,
    CONTACT_CUTOFF_A,
    COORDINATING_ELEMENTS,
    HBOND_DISTANCE_CUTOFF_A,
    HBOND_DONOR_ACCEPTOR_ELEMENTS,
    HBOND_DONOR_HYDROGEN_MAX_DISTANCE_A,
    HBOND_MIN_ANGLE_DEG,
    HYDROPHOBIC_CUTOFF_A,
    HYDROPHOBIC_RESNAMES,
    METAL_COORDINATION_CUTOFF_A,
    METAL_ELEMENTS,
    MOLINSPECT_HEURISTIC_BACKEND,
    NEGATIVE_ATOMS,
    PI_STACKING_CENTER_CUTOFF_A,
    PI_STACKING_MAX_ANGLE_DEVIATION_DEG,
    PI_STACKING_MAX_PARALLEL_OFFSET_A,
    POSITIVE_ATOMS,
    SALT_BRIDGE_CUTOFF_A,
    STERIC_CLASH_MIN_VDW_OVERLAP_A,
    WATER_BRIDGE_MAX_DISTANCE_A,
    definition,
    relation_priority_for_type as _definition_relation_priority_for_type,
)


@dataclass(frozen=True, slots=True)
class RelationSemantics:
    """Classified relation details for the closest atom pair."""

    relation_type: str
    category: str
    confidence: str
    definition_id: str
    definition_source: str
    reference_keys: tuple[str, ...]
    cutoff_A: float
    angle_deg: float | None = None
    angle_metric: str | None = None
    extra_evidence: tuple[EvidenceItem, ...] = ()
    limitations: tuple[str, ...] = ()


def relation_priority_for_type(relation_type: str) -> int:
    """Return stable sort priority for public relation labels."""

    return _definition_relation_priority_for_type(relation_type)


def relation_for_atomgroups(
    source_atoms: Any,
    target_atoms: Any,
    source_id: str,
    target_id: str,
    pair: Any,
) -> Relation:
    """Build a public relation model from closest-pair geometry."""

    semantics = classify_relation(source_atoms, target_atoms, pair)
    source_atom = atom_label(pair.atom_a)
    target_atom = atom_label(pair.atom_b)
    evidence = [
        EvidenceItem(
            type="metric",
            metric="closest_heavy_atom_distance",
            value=pair.distance_A,
            unit="angstrom",
            source=f"{source_atom}->{target_atom}",
        )
    ]
    if semantics.angle_deg is not None:
        evidence.append(
            EvidenceItem(
                type="metric",
                metric=semantics.angle_metric or "angle",
                value=semantics.angle_deg,
                unit="degree",
                source=f"{source_atom}->{target_atom}",
            )
        )
    evidence.extend(semantics.extra_evidence)

    return Relation(
        source=source_id,
        target=target_id,
        type=semantics.relation_type,
        category=semantics.category,
        confidence=semantics.confidence,
        backend=definition(semantics.definition_id).backend or MOLINSPECT_HEURISTIC_BACKEND,
        definition_id=semantics.definition_id,
        definition_source=semantics.definition_source,
        reference_keys=list(semantics.reference_keys),
        min_distance_A=pair.distance_A,
        cutoff_A=semantics.cutoff_A,
        angle_deg=semantics.angle_deg,
        source_atom=source_atom,
        target_atom=target_atom,
        evidence=evidence,
        limitations=list(semantics.limitations),
    )


def water_bridge_relation(
    source_id: str,
    target_id: str,
    direct_pair: Any,
    source_water_pair: Any,
    target_water_pair: Any,
    water_id: str,
) -> Relation:
    """Build a compact relation for one source-water-target polar bridge."""

    source_atom = atom_label(source_water_pair.atom_a)
    water_atom_from_source = atom_label(source_water_pair.atom_b)
    target_atom = atom_label(target_water_pair.atom_a)
    water_atom_from_target = atom_label(target_water_pair.atom_b)
    return Relation(
        source=source_id,
        target=target_id,
        type="water_bridge_candidate",
        category="water_mediated",
        confidence="candidate",
        backend=definition("water_bridge_candidate").backend or MOLINSPECT_HEURISTIC_BACKEND,
        definition_id="water_bridge_candidate",
        definition_source=definition("water_bridge_candidate").source,
        reference_keys=list(definition("water_bridge_candidate").reference_keys),
        min_distance_A=direct_pair.distance_A,
        cutoff_A=WATER_BRIDGE_MAX_DISTANCE_A,
        source_atom=source_atom,
        target_atom=target_atom,
        evidence=[
            EvidenceItem(
                type="metric",
                metric="direct_closest_heavy_atom_distance",
                value=direct_pair.distance_A,
                unit="angstrom",
                source=f"{source_atom}->{target_atom}",
            ),
            EvidenceItem(
                type="metric",
                metric="source_to_water_distance",
                value=source_water_pair.distance_A,
                unit="angstrom",
                source=f"{source_atom}->{water_atom_from_source}",
            ),
            EvidenceItem(
                type="metric",
                metric="water_to_target_distance",
                value=target_water_pair.distance_A,
                unit="angstrom",
                source=f"{water_atom_from_target}->{target_atom}",
            ),
            EvidenceItem(
                type="object",
                metric="bridging_water",
                value=water_id,
                source=water_id,
            ),
        ],
        limitations=list(definition("water_bridge_candidate").limitations),
    )


def classify_relation(source_atoms: Any, target_atoms: Any, pair: Any) -> RelationSemantics:
    """Classify closest-pair relation using compact, documented geometric rules."""

    distance_A = pair.distance_A
    if _metal_coordination_pair(pair) and distance_A <= METAL_COORDINATION_CUTOFF_A:
        return RelationSemantics(
            relation_type="metal_coordination",
            category="coordination",
            confidence="geometry",
            definition_id="metal_coordination",
            definition_source=definition("metal_coordination").source,
            reference_keys=definition("metal_coordination").reference_keys,
            cutoff_A=METAL_COORDINATION_CUTOFF_A,
            limitations=definition("metal_coordination").limitations,
        )

    bond_provenance = bond_provenance_for_pair(pair)
    if bond_provenance is not None:
        definition_record = definition(bond_provenance.definition_id)
        return RelationSemantics(
            relation_type=bond_provenance.definition_id,
            category="covalent",
            confidence=bond_provenance.confidence,
            definition_id=bond_provenance.definition_id,
            definition_source=definition_record.source,
            reference_keys=definition_record.reference_keys,
            cutoff_A=bond_provenance.cutoff_A or 0.0,
            extra_evidence=(
                EvidenceItem(
                    type="method",
                    metric="bond_provenance",
                    value=bond_provenance.method,
                    source=bond_provenance.definition_id,
                ),
            ),
            limitations=definition_record.limitations,
        )

    overlap_A = vdw_overlap_A(pair.atom_a, pair.atom_b, distance_A)
    if (
        overlap_A is not None
        and overlap_A >= STERIC_CLASH_MIN_VDW_OVERLAP_A
        and not _pair_has_metal(pair)
    ):
        return RelationSemantics(
            relation_type="steric_clash",
            category="steric",
            confidence="geometry",
            definition_id="steric_clash",
            definition_source=definition("steric_clash").source,
            reference_keys=definition("steric_clash").reference_keys,
            cutoff_A=STERIC_CLASH_MIN_VDW_OVERLAP_A,
            extra_evidence=(
                EvidenceItem(
                    type="metric",
                    metric="vdw_overlap",
                    value=overlap_A,
                    unit="angstrom",
                    source=f"{atom_label(pair.atom_a)}->{atom_label(pair.atom_b)}",
                ),
            ),
            limitations=definition("steric_clash").limitations,
        )
    if _charged_atom_pair(pair) and distance_A <= SALT_BRIDGE_CUTOFF_A:
        return RelationSemantics(
            relation_type="salt_bridge",
            category="electrostatic",
            confidence="geometry",
            definition_id="salt_bridge",
            definition_source=definition("salt_bridge").source,
            reference_keys=definition("salt_bridge").reference_keys,
            cutoff_A=SALT_BRIDGE_CUTOFF_A,
            limitations=definition("salt_bridge").limitations,
        )

    hbond_angle = _hydrogen_bond_angle_deg(source_atoms, target_atoms, pair)
    if hbond_angle is not None and distance_A <= HBOND_DISTANCE_CUTOFF_A:
        return RelationSemantics(
            relation_type="hydrogen_bond",
            category="polar",
            confidence="geometry",
            definition_id="hydrogen_bond",
            definition_source=definition("hydrogen_bond").source,
            reference_keys=definition("hydrogen_bond").reference_keys,
            cutoff_A=HBOND_DISTANCE_CUTOFF_A,
            angle_deg=hbond_angle,
            angle_metric="donor_hydrogen_acceptor_angle",
        )
    if _pair_has_hbond_elements(pair) and distance_A <= HBOND_DISTANCE_CUTOFF_A:
        return RelationSemantics(
            relation_type="polar_contact_candidate",
            category="polar",
            confidence="candidate",
            definition_id="polar_contact_candidate",
            definition_source=definition("polar_contact_candidate").source,
            reference_keys=definition("polar_contact_candidate").reference_keys,
            cutoff_A=HBOND_DISTANCE_CUTOFF_A,
            limitations=definition("polar_contact_candidate").limitations,
        )

    pi_geometry = _pi_stacking_geometry(pair)
    if pi_geometry is not None:
        center_distance_A, angle_deg, parallel_offset_A = pi_geometry
        return RelationSemantics(
            relation_type="pi_stacking",
            category="aromatic",
            confidence="geometry",
            definition_id="pi_stacking",
            definition_source=definition("pi_stacking").source,
            reference_keys=definition("pi_stacking").reference_keys,
            cutoff_A=PI_STACKING_CENTER_CUTOFF_A,
            angle_deg=angle_deg,
            angle_metric="ring_plane_angle",
            extra_evidence=(
                EvidenceItem(
                    type="metric",
                    metric="ring_center_distance",
                    value=center_distance_A,
                    unit="angstrom",
                    source=f"{atom_label(pair.atom_a)}->{atom_label(pair.atom_b)}",
                ),
                EvidenceItem(
                    type="metric",
                    metric="parallel_ring_offset",
                    value=parallel_offset_A,
                    unit="angstrom",
                    source=f"{atom_label(pair.atom_a)}->{atom_label(pair.atom_b)}",
                ),
            ),
            limitations=definition("pi_stacking").limitations,
        )
    if _hydrophobic_pair(pair) and distance_A <= HYDROPHOBIC_CUTOFF_A:
        return RelationSemantics(
            relation_type="hydrophobic_contact",
            category="hydrophobic",
            confidence="geometry",
            definition_id="hydrophobic_contact",
            definition_source=definition("hydrophobic_contact").source,
            reference_keys=definition("hydrophobic_contact").reference_keys,
            cutoff_A=HYDROPHOBIC_CUTOFF_A,
        )
    if distance_A <= CONTACT_CUTOFF_A:
        return RelationSemantics(
            relation_type="nonbonded_contact",
            category="generic_contact",
            confidence="geometry",
            definition_id="nonbonded_contact",
            definition_source=definition("nonbonded_contact").source,
            reference_keys=definition("nonbonded_contact").reference_keys,
            cutoff_A=CONTACT_CUTOFF_A,
            limitations=definition("nonbonded_contact").limitations,
        )
    return RelationSemantics(
        relation_type="near",
        category="proximity",
        confidence="geometry",
        definition_id="near",
        definition_source=definition("near").source,
        reference_keys=definition("near").reference_keys,
        cutoff_A=CONTACT_CUTOFF_A,
    )


def atom_label(atom: Any) -> str:
    """Return compact evidence label for an atom."""

    residue = atom.residue
    chain = chain_for_residue(residue)
    return f"{chain}:{residue.resid}:{residue.resname}:{atom.name}"


def chain_for_residue(residue: Any) -> str:
    for atom in residue.atoms:
        chain = str(getattr(atom, "chainID", "")).strip()
        if chain:
            return chain
    segment = str(getattr(residue, "segid", "")).strip()
    return segment or "_"


def _pair_has_metal(pair: Any) -> bool:
    return atom_element(pair.atom_a) in METAL_ELEMENTS or atom_element(pair.atom_b) in METAL_ELEMENTS


def _metal_coordination_pair(pair: Any) -> bool:
    element_a = atom_element(pair.atom_a)
    element_b = atom_element(pair.atom_b)
    return (element_a in METAL_ELEMENTS and element_b in COORDINATING_ELEMENTS) or (
        element_b in METAL_ELEMENTS and element_a in COORDINATING_ELEMENTS
    )


def _pair_has_hbond_elements(pair: Any) -> bool:
    return (
        atom_element(pair.atom_a) in HBOND_DONOR_ACCEPTOR_ELEMENTS
        and atom_element(pair.atom_b) in HBOND_DONOR_ACCEPTOR_ELEMENTS
    )


def _charged_atom_pair(pair: Any) -> bool:
    return (_atom_charge_class(pair.atom_a), _atom_charge_class(pair.atom_b)) in {
        ("positive", "negative"),
        ("negative", "positive"),
    }


def _atom_charge_class(atom: Any) -> str | None:
    residue = atom.residue
    resname = str(getattr(residue, "resname", "")).strip().upper()
    atom_name = str(getattr(atom, "name", "")).strip().upper()
    if atom_name in POSITIVE_ATOMS.get(resname, set()):
        return "positive"
    if atom_name in NEGATIVE_ATOMS.get(resname, set()):
        return "negative"
    return None


def _hydrogen_bond_angle_deg(source_atoms: Any, target_atoms: Any, pair: Any) -> float | None:
    candidates = [
        _angle_for_donor_acceptor(pair.atom_a, pair.atom_b, source_atoms),
        _angle_for_donor_acceptor(pair.atom_b, pair.atom_a, target_atoms),
    ]
    valid_angles = [angle for angle in candidates if angle is not None and angle >= HBOND_MIN_ANGLE_DEG]
    if not valid_angles:
        return None
    return round(max(valid_angles), 1)


def _angle_for_donor_acceptor(donor: Any, acceptor: Any, donor_group: Any) -> float | None:
    if atom_element(donor) not in HBOND_DONOR_ACCEPTOR_ELEMENTS:
        return None
    if atom_element(acceptor) not in HBOND_DONOR_ACCEPTOR_ELEMENTS:
        return None
    donor_position = np.asarray(donor.position, dtype=float)
    acceptor_position = np.asarray(acceptor.position, dtype=float)
    for hydrogen in donor_group:
        if hydrogen.residue is not donor.residue or atom_element(hydrogen) != "H":
            continue
        hydrogen_position = np.asarray(hydrogen.position, dtype=float)
        if np.linalg.norm(hydrogen_position - donor_position) > HBOND_DONOR_HYDROGEN_MAX_DISTANCE_A:
            continue
        angle = _angle_degrees(donor_position, hydrogen_position, acceptor_position)
        if angle is not None:
            return angle
    return None


def _angle_degrees(point_a: np.ndarray, point_b: np.ndarray, point_c: np.ndarray) -> float | None:
    vector_a = point_a - point_b
    vector_c = point_c - point_b
    norm_product = float(np.linalg.norm(vector_a) * np.linalg.norm(vector_c))
    if norm_product == 0.0:
        return None
    cosine = float(np.dot(vector_a, vector_c) / norm_product)
    return float(np.degrees(np.arccos(np.clip(cosine, -1.0, 1.0))))


def _hydrophobic_pair(pair: Any) -> bool:
    return _is_hydrophobic_atom(pair.atom_a) and _is_hydrophobic_atom(pair.atom_b)


def _is_hydrophobic_atom(atom: Any) -> bool:
    resname = str(getattr(atom.residue, "resname", "")).strip().upper()
    return atom_element(atom) == "C" and resname in HYDROPHOBIC_RESNAMES


def _pi_stacking_geometry(pair: Any) -> tuple[float, float, float] | None:
    ring_a = _ring_geometry(pair.atom_a.residue)
    ring_b = _ring_geometry(pair.atom_b.residue)
    if ring_a is None or ring_b is None:
        return None
    center_a, normal_a = ring_a
    center_b, normal_b = ring_b
    center_vector = center_b - center_a
    center_distance_A = float(np.linalg.norm(center_vector))
    if center_distance_A > PI_STACKING_CENTER_CUTOFF_A:
        return None

    normal_dot = float(np.dot(normal_a, normal_b))
    angle_deg = float(np.degrees(np.arccos(np.clip(abs(normal_dot), -1.0, 1.0))))
    is_parallel = angle_deg <= PI_STACKING_MAX_ANGLE_DEVIATION_DEG
    is_t_shaped = abs(angle_deg - 90.0) <= PI_STACKING_MAX_ANGLE_DEVIATION_DEG
    if not (is_parallel or is_t_shaped):
        return None

    normal_separation_A = abs(float(np.dot(center_vector, normal_a)))
    parallel_offset_A = float(max(center_distance_A**2 - normal_separation_A**2, 0.0) ** 0.5)
    if is_parallel and parallel_offset_A > PI_STACKING_MAX_PARALLEL_OFFSET_A:
        return None

    return round(center_distance_A, 3), round(angle_deg, 1), round(parallel_offset_A, 3)


def _ring_geometry(residue: Any) -> tuple[np.ndarray, np.ndarray] | None:
    resname = str(getattr(residue, "resname", "")).strip().upper()
    atom_names = AROMATIC_RING_ATOMS.get(resname)
    if atom_names is None:
        return None
    atoms_by_name = {str(atom.name).strip().upper(): atom for atom in residue.atoms}
    positions = [
        np.asarray(atoms_by_name[name].position, dtype=float)
        for name in atom_names
        if name in atoms_by_name
    ]
    if len(positions) < 3:
        return None
    coordinates = np.asarray(positions, dtype=float)
    center = coordinates.mean(axis=0)
    _, _, vh = np.linalg.svd(coordinates - center)
    normal = vh[-1]
    norm = float(np.linalg.norm(normal))
    if norm == 0.0:
        return None
    return center, normal / norm
