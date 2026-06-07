"""Transparent context-scale presets for common inspection tasks."""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from dataclasses import dataclass
from typing import Any, cast

from ..definitions import DEFAULT_CONTEXT_RADIUS_A, definition
from ..errors import MetricError
from ..notation import BUDGET_NOTATION_HELP, CONTEXT_SCALE_NOTATION_HELP
from ..schemas import (
    ContextBudget,
    ContextFocus,
    ContextScale,
    ContextScaleInfo,
    ContextScalesResult,
)


@dataclass(frozen=True, slots=True)
class ContextScaleSpec:
    """Concrete controls behind one named context scale."""

    id: ContextScale
    radius_A: float
    budget: ContextBudget
    focus: tuple[ContextFocus, ...]
    description: str


@dataclass(frozen=True, slots=True)
class ContextSettings:
    """Effective context controls after applying a scale and explicit overrides."""

    scale: ContextScale | None
    radius_A: float
    budget: ContextBudget
    focus: ContextFocus | Sequence[ContextFocus] | None


CONTEXT_SCALE_DESCRIPTIONS: dict[ContextScale, str] = {
    "chemical_contacts": "Close non-water objects and contact-type relations around the selection.",
    "residue_environment": "Broader local neighborhood around a residue, ligand, site, or region.",
    "ligand_binding_site": "Ligand/ion contact shell and nearby binding-site residues.",
    "metal_coordination": "Metal-coordination-centered local context.",
    "protein_interface": "Inter-chain contacts around a protein chain or interface selection.",
    "broad_environment": "Wide static neighborhood for orientation before narrowing the inspection.",
}


def context_scale_specs() -> tuple[ContextScaleSpec, ...]:
    """Return the supported context scales and their concrete settings."""

    return tuple(CONTEXT_SCALE_SPECS.values())


def context_scales() -> ContextScalesResult:
    """Return user-facing context scale presets and their concrete controls."""

    scales = [
        ContextScaleInfo(
            id=spec.id,
            radius_A=spec.radius_A,
            budget=spec.budget,
            focus=list(spec.focus),
            description=spec.description,
        )
        for spec in context_scale_specs()
    ]
    return ContextScalesResult(
        scales=scales,
        count=len(scales),
        limitations=list(definition("context_scale").limitations),
    )


def resolve_context_settings(
    scale: ContextScale | str | None,
    radius: float | None,
    budget: ContextBudget | str | None,
    focus: ContextFocus | Sequence[ContextFocus] | None,
) -> ContextSettings:
    """Resolve a context scale and explicit overrides into concrete controls."""

    spec = _context_scale_spec(scale)
    effective_radius = (
        radius if radius is not None else spec.radius_A if spec else DEFAULT_CONTEXT_RADIUS_A
    )
    effective_budget = _context_budget(
        budget if budget is not None else spec.budget if spec else "small"
    )
    effective_focus = focus if focus is not None else spec.focus if spec else None
    return ContextSettings(
        scale=spec.id if spec else None,
        radius_A=effective_radius,
        budget=effective_budget,
        focus=effective_focus,
    )


def _context_scale_spec(scale: ContextScale | str | None) -> ContextScaleSpec | None:
    if scale is None:
        return None
    if not isinstance(scale, str):
        raise MetricError(CONTEXT_SCALE_NOTATION_HELP)
    key = scale.strip().lower()
    try:
        return CONTEXT_SCALE_SPECS[cast(ContextScale, key)]
    except KeyError as exc:
        raise MetricError(CONTEXT_SCALE_NOTATION_HELP) from exc


def _context_budget(budget: ContextBudget | str) -> ContextBudget:
    if budget in {"small", "medium", "large"}:
        return cast(ContextBudget, budget)
    raise MetricError(BUDGET_NOTATION_HELP)


def _build_context_scale_specs() -> dict[ContextScale, ContextScaleSpec]:
    specs: dict[ContextScale, ContextScaleSpec] = {}
    for raw_id, raw_spec in definition("context_scale").parameters.items():
        if not isinstance(raw_spec, Mapping):
            raise RuntimeError("context_scale definition parameters must be mappings")
        scale_id = cast(ContextScale, str(raw_id))
        specs[scale_id] = ContextScaleSpec(
            id=scale_id,
            radius_A=float(_required_context_scale_field(raw_spec, "radius_A")),
            budget=_context_budget(str(_required_context_scale_field(raw_spec, "budget"))),
            focus=_context_focus_tuple(_required_context_scale_field(raw_spec, "focus")),
            description=CONTEXT_SCALE_DESCRIPTIONS[scale_id],
        )
    return specs


def _required_context_scale_field(raw_spec: Mapping[Any, Any], field: str) -> Any:
    try:
        return raw_spec[field]
    except KeyError as exc:
        raise RuntimeError(f"context_scale definition missing {field!r}") from exc


def _context_focus_tuple(raw_focus: Any) -> tuple[ContextFocus, ...]:
    if isinstance(raw_focus, str):
        return (cast(ContextFocus, raw_focus),)
    if isinstance(raw_focus, Sequence):
        return tuple(cast(ContextFocus, str(value)) for value in raw_focus)
    raise RuntimeError("context_scale focus must be a string or sequence of strings")


CONTEXT_SCALE_SPECS: dict[ContextScale, ContextScaleSpec] = _build_context_scale_specs()
