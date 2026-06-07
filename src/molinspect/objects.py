"""Stable molecular object extraction."""

from __future__ import annotations

from collections.abc import Iterable, Sequence
from typing import Any

from .errors import ObjectQueryError
from .notation import CONTAINS_NOTATION_HELP, OBJECT_TYPE_NOTATION_HELP
from .schemas import ObjectRef, ObjectsResult, ObjectType
from .definitions import (
    ION_RESNAMES,
    MAX_RESOLVED_OBJECT_REFS,
    WATER_RESNAMES,
    definition,
)

TYPE_ALIASES = {
    "atoms": "atom",
    "residues": "residue",
    "ligands": "ligand",
    "ions": "ion",
    "waters": "water",
    "water": "water",
    "chains": "chain",
    "segments": "chain",
    "segment": "chain",
    "chain/segment": "chain",
    "secondary_structures": "secondary_structure",
    "sse": "secondary_structure",
    "sses": "secondary_structure",
    "loops": "loop",
    "ligand_contact_shells": "ligand_contact_shell",
    "contact_shells": "ligand_contact_shell",
    "pockets": "pocket",
    "interchain_contact_interfaces": "interchain_contact_interface",
    "contact_interfaces": "interchain_contact_interface",
    "biological_interfaces": "biological_interface",
    "pisa_interfaces": "biological_interface",
}
BASE_OBJECT_TYPES = frozenset({"atom", "residue", "ligand", "ion", "water", "chain"})
LANDMARK_OBJECT_TYPES = frozenset(
    {
        "secondary_structure",
        "loop",
        "ligand_contact_shell",
        "pocket",
        "interchain_contact_interface",
        "biological_interface",
    }
)
SUPPORTED_OBJECT_TYPES = BASE_OBJECT_TYPES | LANDMARK_OBJECT_TYPES


def list_objects(
    universe: Any,
    object_type: str | Sequence[str] | None = None,
    contains: str | None = None,
    limit: int = 50,
    extra_objects: Sequence[ObjectRef] = (),
) -> ObjectsResult:
    """Return compact, stable object references for the loaded topology."""

    if limit < 1:
        raise ObjectQueryError("limit must be >= 1")
    literal = _normalize_contains(contains)

    requested_types = _normalize_requested_types(object_type)
    default_without_atoms = object_type is None
    objects = list(_iter_objects(universe, requested_types))
    objects.extend(obj for obj in extra_objects if obj.type in requested_types)
    objects = [obj for obj in objects if literal is None or _matches_literal(obj, literal)]

    limitations: list[str] = []
    if default_without_atoms:
        limitations.append("Atom objects are omitted by default; request type='atom' to include them.")

    returned = objects[:limit]
    return ObjectsResult(
        objects=returned,
        count=len(objects),
        returned=len(returned),
        truncated=len(objects) > limit,
        limitations=limitations,
    )


def object_ids_for_atomgroup(atomgroup: Any) -> list[str]:
    """Represent a selected atom group by residue-level objects when possible."""

    return [ref.id for ref in object_refs_for_atomgroup(atomgroup)]


def object_refs_for_atomgroup(atomgroup: Any) -> list[ObjectRef]:
    """Represent a selected atom group by compact residue-level references."""

    residue_ixs = set(getattr(atomgroup.residues, "ix", []))
    if len(residue_ixs) > MAX_RESOLVED_OBJECT_REFS:
        return [_selection_region_ref(atomgroup, len(residue_ixs))]

    classification_sets = _classification_sets(atomgroup.universe)
    refs: list[ObjectRef] = []
    for residue in atomgroup.universe.residues:
        if residue.ix in residue_ixs:
            refs.append(object_ref_for_residue(residue, classification_sets))
    if refs:
        return refs
    return [_object_ref_for_atom(atom) for atom in atomgroup]


def uses_selection_region(atomgroup: Any) -> bool:
    """Return whether an atom group should be compacted as a selection region."""

    residue_ixs = set(getattr(atomgroup.residues, "ix", []))
    return len(residue_ixs) > MAX_RESOLVED_OBJECT_REFS


def object_ref_for_residue(
    residue: Any,
    classification_sets: tuple[set[int], set[int]] | None = None,
) -> ObjectRef:
    """Return the compact object reference for one residue-like object."""

    return _object_ref_for_residue(residue, residue_object_type(residue, classification_sets))


def object_id_for_atom(atom: Any) -> str:
    residue = atom.residue
    chain = _id_part(_chain_or_segment_for_residue(residue))
    resid = _residue_number_with_icode(residue)
    resname = _id_part(getattr(residue, "resname", None))
    atom_name = _id_part(getattr(atom, "name", None))
    return f"atom:{chain}:{resid}:{resname}:{atom_name}:{atom.ix}"


def object_id_for_residue(
    residue: Any,
    classification_sets: tuple[set[int], set[int]] | None = None,
    object_type: ObjectType | None = None,
) -> str:
    object_type = object_type or residue_object_type(residue, classification_sets)
    chain = _id_part(_chain_or_segment_for_residue(residue))
    resid = residue_number_with_icode(residue)
    resname = _id_part(getattr(residue, "resname", None))
    return f"{object_type}:{chain}:{resid}:{resname}"


def residue_object_type(
    residue: Any,
    classification_sets: tuple[set[int], set[int]] | None = None,
) -> ObjectType:
    resname = str(getattr(residue, "resname", "")).strip().upper()
    if resname in WATER_RESNAMES:
        return "water"
    if resname in ION_RESNAMES and len(residue.atoms) <= definition("ion").parameters["max_atom_count"]:
        return "ion"

    protein_ixs, nucleic_ixs = classification_sets or _classification_sets(residue.universe)
    if residue.ix in protein_ixs or residue.ix in nucleic_ixs:
        return "residue"
    return "ligand"


def chain_identifier_for_atom(atom: Any) -> str:
    chain_id = _clean(getattr(atom, "chainID", None))
    if chain_id:
        return chain_id
    segid = _clean(getattr(atom, "segid", None))
    return segid or "_"


def chain_identifier_for_residue(residue: Any) -> str:
    return _chain_or_segment_for_residue(residue)


def residue_number_with_icode(residue: Any) -> str:
    """Return PDB-style residue number including insertion code when present."""

    resid = str(getattr(residue, "resid", "")).strip() or "_"
    icode = _residue_icode(residue)
    return f"{resid}{icode}" if icode else resid


def water_selection_expression() -> str:
    return _resname_expression(WATER_RESNAMES)


def ion_selection_expression() -> str:
    return _resname_expression(ION_RESNAMES)


def _iter_objects(universe: Any, requested_types: set[str]) -> Iterable[ObjectRef]:
    classification_sets = _classification_sets(universe)

    if "chain" in requested_types:
        yield from _chain_objects(universe)

    residue_types = {"residue", "ligand", "ion", "water"} & requested_types
    if residue_types:
        for residue in universe.residues:
            obj_type = residue_object_type(residue, classification_sets)
            if obj_type in residue_types:
                yield _object_ref_for_residue(residue, obj_type)

    if "atom" in requested_types:
        for atom in universe.atoms:
            yield _object_ref_for_atom(atom)


def _chain_objects(universe: Any) -> Iterable[ObjectRef]:
    atom_counts: dict[str, int] = {}
    segments: dict[str, str] = {}
    for atom in universe.atoms:
        chain = chain_identifier_for_atom(atom)
        atom_counts[chain] = atom_counts.get(chain, 0) + 1
        segments.setdefault(chain, _clean(getattr(atom, "segid", None)))

    for chain, atom_count in atom_counts.items():
        yield ObjectRef(
            id=f"chain:{_id_part(chain)}",
            type="chain",
            name=chain,
            chain=chain,
            segment=segments[chain],
            atom_count=atom_count,
            annotations=_object_definition_annotations("chain"),
        )


def _object_ref_for_residue(residue: Any, obj_type: ObjectType) -> ObjectRef:
    chain = chain_identifier_for_residue(residue)
    segment = _clean(getattr(residue, "segid", None))
    resname = _clean(getattr(residue, "resname", None))
    resid = getattr(residue, "resid", None)
    icode = _residue_icode(residue) or None
    return ObjectRef(
        id=object_id_for_residue(residue, object_type=obj_type),
        type=obj_type,
        name=f"{resname}{resid}" if resname and resid is not None else resname,
        chain=chain,
        segment=segment,
        resid=resid,
        icode=icode,
        resname=resname,
        atom_count=len(residue.atoms),
        annotations=_object_definition_annotations(obj_type),
    )


def _object_ref_for_atom(atom: Any) -> ObjectRef:
    residue = atom.residue
    chain = chain_identifier_for_atom(atom)
    segment = _clean(getattr(atom, "segid", None))
    resname = _clean(getattr(residue, "resname", None))
    resid = getattr(residue, "resid", None)
    icode = _residue_icode(residue) or None
    atom_name = _clean(getattr(atom, "name", None))
    altloc = _atom_altloc(atom)
    return ObjectRef(
        id=object_id_for_atom(atom),
        type="atom",
        name=atom_name,
        chain=chain,
        segment=segment,
        resid=resid,
        icode=icode,
        resname=resname,
        atom_name=atom_name,
        altloc=altloc,
        atom_index=atom.ix,
        atom_count=1,
        annotations=_object_definition_annotations("atom"),
    )


def _normalize_requested_types(object_type: str | Sequence[str] | None) -> set[str]:
    if object_type is None:
        return {"chain", "residue", "ligand", "ion", "water"}

    raw_types: Sequence[str]
    if isinstance(object_type, str):
        raw_types = [object_type]
    else:
        raw_types = object_type

    normalized: set[str] = set()
    for raw_type in raw_types:
        if not isinstance(raw_type, str):
            raise ObjectQueryError(OBJECT_TYPE_NOTATION_HELP)
        value = raw_type.strip().lower()
        value = TYPE_ALIASES.get(value, value)
        if value not in SUPPORTED_OBJECT_TYPES:
            raise ObjectQueryError(f"Unsupported object type: {raw_type!r}. {OBJECT_TYPE_NOTATION_HELP}")
        normalized.add(value)
    return normalized


def _normalize_contains(contains: str | None) -> str | None:
    if contains is None:
        return None
    if not isinstance(contains, str):
        raise ObjectQueryError(CONTAINS_NOTATION_HELP)
    literal = contains.strip()
    if not literal:
        raise ObjectQueryError(CONTAINS_NOTATION_HELP)
    return literal.casefold()


def _matches_literal(obj: ObjectRef, needle: str) -> bool:
    fields = [
        obj.id,
        obj.name,
        obj.chain,
        obj.segment,
        str(obj.resid) if obj.resid is not None else None,
        obj.resname,
        obj.atom_name,
    ]
    return any(field is not None and needle in field.casefold() for field in fields)


def _classification_sets(universe: Any) -> tuple[set[int], set[int]]:
    return _selected_residue_ixs(universe, "protein"), _selected_residue_ixs(universe, "nucleic")


def _selected_residue_ixs(universe: Any, selection: str) -> set[int]:
    try:
        return set(universe.select_atoms(selection).residues.ix)
    except Exception:
        return set()


def _chain_or_segment_for_residue(residue: Any) -> str:
    for atom in residue.atoms:
        chain = _clean(getattr(atom, "chainID", None))
        if chain:
            return chain
    segid = _clean(getattr(residue, "segid", None))
    return segid or "_"


def _residue_number_with_icode(residue: Any) -> str:
    return residue_number_with_icode(residue)


def _residue_icode(residue: Any) -> str:
    icode = _clean(getattr(residue, "icode", None))
    if icode:
        return icode
    try:
        icodes = getattr(residue.atoms, "icodes", [])
    except Exception:
        return ""
    for value in icodes:
        cleaned = _clean(value)
        if cleaned:
            return cleaned
    return ""


def _atom_altloc(atom: Any) -> str | None:
    for attribute in ("altLoc", "altloc"):
        altloc = _clean(getattr(atom, attribute, None))
        if altloc:
            return altloc
    return None


def _resname_expression(resnames: Iterable[str]) -> str:
    return "resname " + " ".join(sorted(resnames))


def _clean(value: Any) -> str:
    if value is None:
        return ""
    cleaned = str(value).strip()
    return "" if cleaned in {"", "None"} else cleaned


def _id_part(value: Any) -> str:
    cleaned = _clean(value) or "_"
    return cleaned.replace(":", "_").replace(" ", "_")


def _selection_region_ref(atomgroup: Any, n_residues: int) -> ObjectRef:
    atom_indices = [int(atom.ix) for atom in atomgroup]
    first_index = min(atom_indices) if atom_indices else 0
    last_index = max(atom_indices) if atom_indices else 0
    fingerprint = _selection_fingerprint(atom_indices)
    return ObjectRef(
        id=f"selection_region:{n_residues}res:{first_index}-{last_index}:{fingerprint}",
        type="selection_region",
        name=f"selection_region:{n_residues}_residues",
        atom_count=len(atomgroup),
        annotations={
            "definition_id": "selection_region",
            "definition_source": definition("selection_region").source,
            "limitations": list(definition("selection_region").limitations),
        },
    )


def _object_definition_annotations(object_type: ObjectType) -> dict[str, Any]:
    try:
        definition_record = definition(object_type)
    except KeyError:
        return {}
    output: dict[str, Any] = {
        "definition_id": definition_record.id,
        "definition_source": definition_record.source,
    }
    if definition_record.reference_keys:
        output["reference_keys"] = list(definition_record.reference_keys)
    if definition_record.limitations:
        output["limitations"] = list(definition_record.limitations)
    return output


def _selection_fingerprint(atom_indices: list[int]) -> str:
    value = 0
    for index in atom_indices:
        value = ((value * 1315423911) ^ index) & 0xFFFFFFFF
    return f"{value:08x}"
