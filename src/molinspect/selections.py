"""Selection translation and resolution."""

from __future__ import annotations

import re
from dataclasses import dataclass
from typing import Any

from .errors import FrameError, SelectionResolutionError
from .notation import FRAME_NOTATION_HELP, SELECTION_NOTATION_HELP
from .objects import (
    ion_selection_expression,
    object_id_for_atom,
    object_ids_for_atomgroup,
    uses_selection_region,
    water_selection_expression,
)
from .schemas import SelectionResult
from .definitions import MAX_RESOLVED_OBJECT_REFS

_CHAIN_PATTERN = re.compile(r"\bchain\s+([A-Za-z0-9_.-]+)\b")
_SHORT_RESIDUE_PATTERN = re.compile(r"^\s*([A-Za-z0-9_.-]+):([+-]?\d+[A-Za-z]?)\s*$")
_AROUND_OF_PATTERN = re.compile(r"^\s*around\s+([0-9]+(?:\.[0-9]+)?)\s+of\s+(.+?)\s*$", re.I)
_RESIDUE_OBJECT_TYPES = frozenset({"residue", "ligand", "ion", "water"})


@dataclass(frozen=True)
class ResolvedSelection:
    """Backend-resolved representation of a public MolInspect selection."""

    selection: str
    expression: str
    frame: int
    atomgroup: Any
    resolved_objects: tuple[str, ...]
    limitations: tuple[str, ...] = ()

    def to_model(self) -> SelectionResult:
        return SelectionResult(
            selection=self.selection,
            expression=self.expression,
            frame=self.frame,
            n_atoms=len(self.atomgroup),
            resolved_objects=list(self.resolved_objects),
            limitations=list(self.limitations),
        )


def normalize_frame(frame: int | str, n_frames: int) -> int:
    """Normalize public frame references to a zero-based frame index."""

    if n_frames < 1:
        raise FrameError("loaded target has no frames")

    if isinstance(frame, str):
        value = frame.strip().lower()
        if value == "first":
            return 0
        if value == "last":
            return n_frames - 1
        if value == "representative":
            raise FrameError(
                "frame='representative' is reserved for later representative-state selection. "
                f"{FRAME_NOTATION_HELP}"
            )
        try:
            frame = int(value)
        except ValueError as exc:
            raise FrameError(FRAME_NOTATION_HELP) from exc

    if isinstance(frame, bool) or not isinstance(frame, int):
        raise FrameError(FRAME_NOTATION_HELP)

    normalized = frame + n_frames if frame < 0 else frame
    if normalized < 0 or normalized >= n_frames:
        raise FrameError(
            f"frame {frame!r} is outside the valid range 0..{n_frames - 1}. "
            f"{FRAME_NOTATION_HELP}"
        )
    return normalized


def resolve_selection(universe: Any, selection: str, frame: int | str = 0) -> ResolvedSelection:
    """Resolve a public MolInspect selection to an MDAnalysis AtomGroup."""

    if not isinstance(selection, str) or not selection.strip():
        raise SelectionResolutionError(
            f"selection must be a non-empty string. {SELECTION_NOTATION_HELP}"
        )

    frame_index = normalize_frame(frame, len(universe.trajectory))
    universe.trajectory[frame_index]

    expression, limitations = translate_selection(universe, selection)
    try:
        atomgroup = universe.select_atoms(expression)
    except Exception as exc:
        raise SelectionResolutionError(
            f"Could not resolve selection {selection!r} as {expression!r}: {exc}. "
            f"{SELECTION_NOTATION_HELP}"
        ) from exc
    if len(atomgroup) == 0:
        raise SelectionResolutionError(
            f"Selection {selection!r} resolved to zero atoms as {expression!r}. "
            f"{SELECTION_NOTATION_HELP}"
        )

    resolved_objects = _resolved_object_ids_for_selection(selection.strip(), atomgroup)
    if uses_selection_region(atomgroup):
        limitations.append(
            "Selection spans more than "
            f"{MAX_RESOLVED_OBJECT_REFS} residue objects; represented as one selection_region."
        )

    return ResolvedSelection(
        selection=selection,
        expression=expression,
        frame=frame_index,
        atomgroup=atomgroup,
        resolved_objects=resolved_objects,
        limitations=tuple(limitations),
    )


def _resolved_object_ids_for_selection(selection: str, atomgroup: Any) -> tuple[str, ...]:
    if selection.startswith("atom:") and len(atomgroup) == 1:
        return (object_id_for_atom(atomgroup[0]),)
    return tuple(object_ids_for_atomgroup(atomgroup))


def translate_selection(universe: Any, selection: str) -> tuple[str, list[str]]:
    """Translate MolInspect-facing selection conveniences to MDAnalysis syntax."""

    expression = _expand_object_id_alias(selection.strip())
    expression = _expand_around_of_alias(_expand_short_residue_alias(expression))
    limitations: list[str] = []

    expression = _CHAIN_PATTERN.sub(lambda match: _chain_replacement(universe, match), expression)

    water_expr = f"({water_selection_expression()})"
    ion_expr = f"({ion_selection_expression()})"
    ligand_expr = f"(not protein and not nucleic and not {water_expr} and not {ion_expr})"

    expression = re.sub(r"\bligand\b", ligand_expr, expression, flags=re.IGNORECASE)
    expression = re.sub(r"\bwater\b", water_expr, expression, flags=re.IGNORECASE)
    expression = re.sub(r"\bion\b", ion_expr, expression, flags=re.IGNORECASE)

    if expression != selection.strip():
        limitations.append("Selection was translated to MDAnalysis syntax.")
    return expression, limitations


def _expand_object_id_alias(selection: str) -> str:
    parts = selection.split(":")
    if len(parts) == 2 and parts[0] == "chain":
        return f"chain {parts[1]}"
    if len(parts) == 6 and parts[0] == "atom":
        atom_index = parts[5]
        if atom_index.isdigit():
            return f"index {atom_index}"
        return selection
    if len(parts) == 4 and parts[0] in _RESIDUE_OBJECT_TYPES:
        _, chain, resid, resname = parts
        selectors: list[str] = []
        if chain != "_":
            selectors.append(f"chain {chain}")
        if resid != "_":
            selectors.append(f"resid {resid}")
        if resname != "_":
            selectors.append(f"resname {resname}")
        return " and ".join(selectors) if selectors else selection
    return selection


def _expand_short_residue_alias(selection: str) -> str:
    match = _SHORT_RESIDUE_PATTERN.match(selection)
    if not match:
        return selection
    chain, resid = match.groups()
    return f"chain {chain} and resid {resid}"


def _expand_around_of_alias(selection: str) -> str:
    match = _AROUND_OF_PATTERN.match(selection)
    if not match:
        return selection
    radius, inner = match.groups()
    inner = inner.strip()
    if inner.startswith("(") and inner.endswith(")"):
        return f"around {radius} {inner}"
    return f"around {radius} ({inner})"


def _chain_replacement(universe: Any, match: re.Match[str]) -> str:
    requested = match.group(1)
    if _has_matching_chain_id(universe, requested):
        return f"chainID {requested}"
    if _has_matching_segment_id(universe, requested):
        return f"segid {requested}"
    return f"chainID {requested}"


def _has_matching_chain_id(universe: Any, requested: str) -> bool:
    try:
        return requested in {
            str(value).strip()
            for value in universe.atoms.chainIDs
            if str(value).strip()
        }
    except Exception:
        return False


def _has_matching_segment_id(universe: Any, requested: str) -> bool:
    try:
        return requested in {str(value).strip() for value in universe.atoms.segids}
    except Exception:
        return False
