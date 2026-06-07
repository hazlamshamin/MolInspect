"""Biological landmark objects derived from structural annotations."""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

from .annotations import AnnotationStore
from .backends.interfaces import (
    detect_biological_interfaces,
    object_ref_for_biological_interface,
    selection_expression_for_biological_interface,
)
from .objects import chain_identifier_for_residue, object_id_for_residue, residue_object_type
from .backends.pockets import (
    detect_pockets,
    object_ref_for_pocket,
    selection_expression_for_pocket,
)
from .schemas import ObjectRef
from .definitions import INTERFACE_DISTANCE_A, LIGAND_CONTACT_DISTANCE_A, definition

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


@dataclass(slots=True)
class LandmarkIndex:
    """Computed landmark objects and residue memberships for one topology."""

    objects: list[ObjectRef]
    selection_by_id: dict[str, str]
    memberships_by_residue_ix: dict[int, set[str]] = field(default_factory=dict)
    secondary_element_by_residue_ix: dict[int, str] = field(default_factory=dict)
    ligand_contact_shells_by_residue_ix: dict[int, set[str]] = field(default_factory=dict)
    pockets_by_residue_ix: dict[int, set[str]] = field(default_factory=dict)
    interchain_contact_interfaces_by_residue_ix: dict[int, set[str]] = field(default_factory=dict)
    biological_interfaces_by_residue_ix: dict[int, set[str]] = field(default_factory=dict)
    limitations: tuple[str, ...] = ()

    def profile_update_for_residue(self, residue_ix: int) -> dict[str, Any]:
        """Return structural-profile fields for one residue's landmark memberships."""

        memberships = sorted(self.memberships_by_residue_ix.get(residue_ix, set()))
        return {
            "landmark_memberships": memberships,
            "secondary_structure_element": self.secondary_element_by_residue_ix.get(residue_ix),
            "ligand_contact_shell_ids": sorted(
                self.ligand_contact_shells_by_residue_ix.get(residue_ix, set())
            ),
            "pocket_ids": sorted(self.pockets_by_residue_ix.get(residue_ix, set())),
            "interchain_contact_interface_ids": sorted(
                self.interchain_contact_interfaces_by_residue_ix.get(residue_ix, set())
            ),
            "biological_interface_ids": sorted(
                self.biological_interfaces_by_residue_ix.get(residue_ix, set())
            ),
        }


def build_landmark_index(
    universe: Any,
    annotations: AnnotationStore,
    include_backend_pockets: bool = False,
    include_backend_interfaces: bool = False,
    source_files: tuple[Path, ...] = (),
) -> LandmarkIndex:
    """Build high-level biological objects from residue annotations."""

    index = LandmarkIndex(objects=[], selection_by_id={})
    _add_secondary_structure_landmarks(universe, annotations, index)
    _add_ligand_contact_shell_landmarks(universe, annotations, index)
    _add_interface_landmarks(universe, annotations, index)
    if include_backend_pockets:
        _add_backend_pocket_landmarks(universe, source_files, index)
    if include_backend_interfaces:
        _add_backend_biological_interface_landmarks(universe, source_files, index)
    index.objects.sort(key=_landmark_sort_key)
    return index


def _add_secondary_structure_landmarks(
    universe: Any,
    annotations: AnnotationStore,
    index: LandmarkIndex,
) -> None:
    run: list[Any] = []
    run_kind: str | None = None
    run_chain: str | None = None

    def flush() -> None:
        nonlocal run, run_kind, run_chain
        if run and run_kind and run_chain:
            _record_secondary_run(run, run_kind, run_chain, annotations, index)
        run = []
        run_kind = None
        run_chain = None

    for residue in universe.residues:
        if residue_object_type(residue) != "residue":
            flush()
            continue
        annotation = annotations.residue(residue)
        kind = annotation.secondary_structure
        chain = chain_identifier_for_residue(residue)
        if not kind:
            flush()
            continue
        if run and (kind != run_kind or chain != run_chain or not _is_next_residue(run[-1], residue)):
            flush()
        run.append(residue)
        run_kind = kind
        run_chain = chain
    flush()


def _record_secondary_run(
    residues: list[Any],
    kind: str,
    chain: str,
    annotations: AnnotationStore,
    index: LandmarkIndex,
) -> None:
    start = _residue_number(residues[0])
    end = _residue_number(residues[-1])
    object_type = "loop" if kind == "loop" else "secondary_structure"
    kind_part = "" if object_type == "loop" else f":{_id_part(kind)}"
    object_id = f"{object_type}:{_id_part(chain)}{kind_part}:{_id_part(start)}-{_id_part(end)}"
    residue_ids = [object_id_for_residue(residue) for residue in residues]
    source = _first_secondary_source(residues, annotations)
    obj = ObjectRef(
        id=object_id,
        type=object_type,  # type: ignore[arg-type]
        name=f"{kind} {chain}:{start}-{end}",
        chain=chain,
        resid=f"{start}-{end}",
        atom_count=sum(len(residue.atoms) for residue in residues),
        annotations={
            "kind": kind,
            "source": source,
            "definition_id": object_type,
            "definition_source": definition(object_type).source,
            "reference_keys": list(definition(object_type).reference_keys),
            "residue_ids": residue_ids,
            "member_count": len(residue_ids),
        },
    )
    index.objects.append(obj)
    index.selection_by_id[object_id] = _or_expression(_residue_expression(residue) for residue in residues)
    for residue in residues:
        index.memberships_by_residue_ix.setdefault(residue.ix, set()).add(object_id)
        index.secondary_element_by_residue_ix[residue.ix] = object_id


def _add_ligand_contact_shell_landmarks(
    universe: Any,
    annotations: AnnotationStore,
    index: LandmarkIndex,
) -> None:
    ligand_residue_map = _residues_by_object_id(universe)
    lining_by_ligand: dict[str, list[Any]] = {}
    for residue in universe.residues:
        annotation = annotations.residue(residue)
        for ligand_id in annotation.ligand_contact_ids:
            lining_by_ligand.setdefault(ligand_id, []).append(residue)

    for ligand_id, lining_residues in sorted(lining_by_ligand.items()):
        ligand_residue = ligand_residue_map.get(ligand_id)
        ligand_expression = _residue_expression(ligand_residue) if ligand_residue is not None else None
        lining_ids = [object_id_for_residue(residue) for residue in lining_residues]
        ligand_parts = _object_id_tail(ligand_id)
        ligand_chain, ligand_resid, ligand_resname = _parse_residue_object_tail(ligand_parts)
        shell_id = f"ligand_contact_shell:{ligand_parts}"
        shell_expression = _or_expression(
            [expr for expr in [ligand_expression, *(_residue_expression(residue) for residue in lining_residues)] if expr]
        )
        shell_definition = definition("ligand_contact_shell")
        shell_annotations = {
            "definition_id": "ligand_contact_shell",
            "definition_source": shell_definition.source,
            "reference_keys": list(shell_definition.reference_keys),
            "ligand_id": ligand_id,
            "lining_residues": lining_ids,
            "lining_residue_count": len(lining_ids),
            "contact_cutoff_A": LIGAND_CONTACT_DISTANCE_A,
            "method": shell_definition.parameters["method"],
            "limitation": shell_definition.limitations[0],
        }
        index.objects.append(
            ObjectRef(
                id=shell_id,
                type="ligand_contact_shell",
                name=f"contact_shell_for_{ligand_resname}{ligand_resid}",
                chain=ligand_chain,
                resid=ligand_resid,
                resname=ligand_resname,
                annotations=shell_annotations,
            )
        )
        index.selection_by_id[shell_id] = shell_expression
        for residue in lining_residues:
            index.memberships_by_residue_ix.setdefault(residue.ix, set()).add(shell_id)
            index.ligand_contact_shells_by_residue_ix.setdefault(residue.ix, set()).add(shell_id)


def _add_interface_landmarks(
    universe: Any,
    annotations: AnnotationStore,
    index: LandmarkIndex,
) -> None:
    residues_by_pair: dict[tuple[str, str], list[Any]] = {}
    for residue in universe.residues:
        own_chain = chain_identifier_for_residue(residue)
        annotation = annotations.residue(residue)
        for partner_chain in annotation.interface_chains:
            chain_a, chain_b = sorted((own_chain, partner_chain))
            pair = (chain_a, chain_b)
            residues_by_pair.setdefault(pair, []).append(residue)

    for chains, residues in sorted(residues_by_pair.items()):
        chain_a, chain_b = chains
        object_id = f"interchain_contact_interface:{_id_part(chain_a)}-{_id_part(chain_b)}"
        residue_ids = sorted({object_id_for_residue(residue) for residue in residues})
        interface_definition = definition("interchain_contact_interface")
        index.objects.append(
            ObjectRef(
                id=object_id,
                type="interchain_contact_interface",
                name=f"interchain_contact_interface {chain_a}-{chain_b}",
                annotations={
                    "definition_id": "interchain_contact_interface",
                    "definition_source": interface_definition.source,
                    "reference_keys": list(interface_definition.reference_keys),
                    "chains": [chain_a, chain_b],
                    "participating_residues": residue_ids,
                    "participating_residue_count": len(residue_ids),
                    "contact_cutoff_A": INTERFACE_DISTANCE_A,
                    "method": interface_definition.parameters["method"],
                    "limitation": interface_definition.limitations[0],
                },
            )
        )
        index.selection_by_id[object_id] = _or_expression(_residue_expression(residue) for residue in residues)
        for residue in residues:
            index.memberships_by_residue_ix.setdefault(residue.ix, set()).add(object_id)
            index.interchain_contact_interfaces_by_residue_ix.setdefault(residue.ix, set()).add(object_id)


def _add_backend_pocket_landmarks(
    universe: Any,
    source_files: tuple[Path, ...],
    index: LandmarkIndex,
) -> None:
    pocket_store = detect_pockets(universe, source_files)
    index.limitations = tuple(dict.fromkeys((*index.limitations, *pocket_store.limitations)))
    residue_lookup = _residues_by_object_id(universe)
    for record in pocket_store.records:
        obj = object_ref_for_pocket(record)
        index.objects.append(obj)
        index.selection_by_id[obj.id] = selection_expression_for_pocket(record)
        for residue_id in record.residue_ids:
            residue = residue_lookup.get(residue_id)
            if residue is None:
                continue
            index.memberships_by_residue_ix.setdefault(residue.ix, set()).add(obj.id)
            index.pockets_by_residue_ix.setdefault(residue.ix, set()).add(obj.id)


def _add_backend_biological_interface_landmarks(
    universe: Any,
    source_files: tuple[Path, ...],
    index: LandmarkIndex,
) -> None:
    interface_store = detect_biological_interfaces(universe, source_files)
    index.limitations = tuple(dict.fromkeys((*index.limitations, *interface_store.limitations)))
    residue_lookup = _residues_by_object_id(universe)
    for record in interface_store.records:
        obj = object_ref_for_biological_interface(record)
        index.objects.append(obj)
        index.selection_by_id[obj.id] = selection_expression_for_biological_interface(record)
        for residue_id in record.residue_ids:
            residue = residue_lookup.get(residue_id)
            if residue is None:
                continue
            index.memberships_by_residue_ix.setdefault(residue.ix, set()).add(obj.id)
            index.biological_interfaces_by_residue_ix.setdefault(residue.ix, set()).add(obj.id)


def _residues_by_object_id(universe: Any) -> dict[str, Any]:
    return {object_id_for_residue(residue): residue for residue in universe.residues}


def _first_secondary_source(residues: list[Any], annotations: AnnotationStore) -> str | None:
    for residue in residues:
        source = annotations.residue(residue).secondary_structure_source
        if source:
            return source
    return None


def _is_next_residue(previous: Any, current: Any) -> bool:
    try:
        return int(current.resid) == int(previous.resid) + 1
    except (TypeError, ValueError):
        return False


def _residue_expression(residue: Any) -> str:
    chain_expr = _chain_expression_for_residue(residue)
    resid = _residue_number(residue)
    resname = str(getattr(residue, "resname", "")).strip()
    parts = [f"resid {resid}"]
    if chain_expr:
        parts.insert(0, chain_expr)
    if resname:
        parts.append(f"resname {resname}")
    return " and ".join(parts)


def _chain_expression_for_residue(residue: Any) -> str:
    for atom in residue.atoms:
        chain = str(getattr(atom, "chainID", "")).strip()
        if chain:
            return f"chainID {chain}"
    segment = str(getattr(residue, "segid", "")).strip()
    return f"segid {segment}" if segment else ""


def _or_expression(expressions: Any) -> str:
    values = [f"({expression})" for expression in expressions if expression]
    return " or ".join(values)


def _residue_number(residue: Any) -> str:
    resid = str(getattr(residue, "resid", "")).strip() or "_"
    icode = str(getattr(residue, "icode", "")).strip()
    return f"{resid}{icode}" if icode else resid


def _object_id_tail(object_id: str) -> str:
    return object_id.split(":", 1)[1] if ":" in object_id else object_id


def _parse_residue_object_tail(tail: str) -> tuple[str | None, str | None, str | None]:
    parts = tail.split(":")
    chain = parts[0] if len(parts) > 0 and parts[0] != "_" else None
    resid = parts[1] if len(parts) > 1 and parts[1] != "_" else None
    resname = parts[2] if len(parts) > 2 and parts[2] != "_" else None
    return chain, resid, resname


def _id_part(value: Any) -> str:
    cleaned = str(value).strip() or "_"
    return cleaned.replace(":", "_").replace(" ", "_")


def _landmark_sort_key(obj: ObjectRef) -> tuple[str, int, str]:
    rank = obj.annotations.get("rank")
    if isinstance(rank, int):
        return (obj.type, rank, obj.id)
    return (obj.type, 10**9, obj.id)
