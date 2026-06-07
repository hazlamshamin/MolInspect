"""Input validation helpers for inspection primitives."""

from __future__ import annotations

from ..definitions import BUDGET_LIMITS
from ..errors import MetricError
from ..notation import BUDGET_NOTATION_HELP


def budget_limit(budget: str) -> int:
    """Return the maximum context objects allowed by a budget label."""

    if not isinstance(budget, str):
        raise MetricError(BUDGET_NOTATION_HELP)
    try:
        return BUDGET_LIMITS[budget]
    except KeyError as exc:
        raise MetricError(BUDGET_NOTATION_HELP) from exc
