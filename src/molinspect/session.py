"""Inspection session facade."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Any, Sequence

from .inspection import compare, context, locate, timeline
from .inspection.scales import context_scales as available_context_scales
from .inspection.timeline_specs import timeline_metrics as available_timeline_metrics
from .objects import list_objects
from .schemas import (
    CompareResult,
    ContextBudget,
    ContextFocus,
    ContextResult,
    ContextScale,
    ContextScalesResult,
    LoadResult,
    LocateResult,
    ObjectsResult,
    SelectionResult,
    TimelineMetricsResult,
    TimelineResult,
)
from .selections import normalize_frame
from .world import InspectionWorld


@dataclass(slots=True)
class InspectionSession:
    """Stateful molecular inspection target.

    The session intentionally exposes a small public surface. Backend details such as
    MDAnalysis remain available internally but should not become LLM-facing API sprawl.
    """

    universe: Any
    world: InspectionWorld
    target_id: str
    name: str | None
    source_files: tuple[Path, ...]
    load_result: LoadResult

    @property
    def n_frames(self) -> int:
        return len(self.universe.trajectory)

    def summary(self) -> LoadResult:
        """Return compact load metadata."""

        return self.load_result

    def normalize_frame(self, frame: int | str = 0) -> int:
        """Normalize a public frame reference for this target."""

        return normalize_frame(frame, self.n_frames)

    def resolve_selection(self, selection: str, frame: int | str = 0) -> SelectionResult:
        """Resolve a public selection and return a compact diagnostic model."""

        return self.world.resolve_selection(selection, frame=frame).to_model()

    def objects(
        self,
        type: str | Sequence[str] | None = None,
        contains: str | None = None,
        frame: int | None = None,
        limit: int = 50,
    ) -> ObjectsResult:
        """List stable molecular objects, optionally filtered by literal text."""

        if frame is not None:
            self.universe.trajectory[self.normalize_frame(frame)]
        include_backend_pockets = _requests_backend_pockets(type)
        include_backend_interfaces = _requests_backend_interfaces(type)
        result = list_objects(
            self.universe,
            object_type=type,
            contains=contains,
            limit=limit,
            extra_objects=self.world.landmark_objects(
                include_backend_pockets=include_backend_pockets,
                include_backend_interfaces=include_backend_interfaces,
            ),
        )
        backend_limitations = self.world.landmark_limitations(
            include_backend_pockets=include_backend_pockets,
            include_backend_interfaces=include_backend_interfaces,
        )
        if not backend_limitations:
            return result
        return result.model_copy(
            update={
                "limitations": list(dict.fromkeys((*result.limitations, *backend_limitations)))
            }
        )

    def locate(
        self,
        selection: str,
        frame: int | str = 0,
        include_metrics: bool = True,
    ) -> LocateResult:
        """Answer where a selected object or region is."""

        return locate(
            self.world,
            selection=selection,
            frame=frame,
            include_metrics=include_metrics,
        )

    def context(
        self,
        selection: str,
        frame: int | str = 0,
        radius: float | None = None,
        budget: ContextBudget | str | None = None,
        focus: ContextFocus | Sequence[ContextFocus] | None = None,
        scale: ContextScale | str | None = None,
    ) -> ContextResult:
        """Return compact context around a selection using explicit controls or a scale."""

        return context(
            self.world,
            selection=selection,
            frame=frame,
            radius=radius,
            budget=budget,
            focus=focus,
            scale=scale,
        )

    def context_scales(self) -> ContextScalesResult:
        """Return supported context-scale presets and their concrete controls."""

        return available_context_scales()

    def timeline(
        self,
        metric: str,
        selection: str | None = None,
        selection1: str | None = None,
        selection2: str | None = None,
        frames: str | tuple[int, int] = "all",
        stride: int | None = None,
    ) -> TimelineResult:
        """Summarize a structural metric, relation, or motion signal over frames."""

        return timeline(
            self.world,
            metric=metric,
            selection=selection,
            selection1=selection1,
            selection2=selection2,
            frames=frames,
            stride=stride,
        )

    def timeline_metrics(self) -> TimelineMetricsResult:
        """Return supported timeline metrics with required arguments and output keys."""

        return available_timeline_metrics()

    def compare(
        self,
        selection: str,
        frame_a: int | str,
        frame_b: int | str,
        radius: float = 8.0,
    ) -> CompareResult:
        """Compare compact local context around a selection between two frames."""

        return compare(
            self.world,
            selection=selection,
            frame_a=frame_a,
            frame_b=frame_b,
            radius=radius,
        )


def _requests_backend_pockets(object_type: str | Sequence[str] | None) -> bool:
    if object_type is None:
        return False
    if isinstance(object_type, str):
        requested = [object_type]
    else:
        requested = list(object_type)
    return any(str(value).strip().lower() in {"pocket", "pockets"} for value in requested)


def _requests_backend_interfaces(object_type: str | Sequence[str] | None) -> bool:
    if object_type is None:
        return False
    if isinstance(object_type, str):
        requested = [object_type]
    else:
        requested = list(object_type)
    return any(
        str(value).strip().lower()
        in {"biological_interface", "biological_interfaces", "pisa_interface", "pisa_interfaces"}
        for value in requested
    )
