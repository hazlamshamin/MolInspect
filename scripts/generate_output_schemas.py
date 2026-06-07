"""Generate the bundled JSON schema artifact from Pydantic models."""

from __future__ import annotations

import json
from pathlib import Path
from typing import Union

from pydantic import TypeAdapter

from molinspect.schemas import (
    CompareResult,
    ContextResult,
    ContextScalesResult,
    LoadResult,
    LocateResult,
    ObjectsResult,
    SelectionResult,
    TimelineMetricsResult,
    TimelineResult,
)

OutputModel = Union[
    LoadResult,
    ObjectsResult,
    SelectionResult,
    LocateResult,
    ContextResult,
    ContextScalesResult,
    TimelineMetricsResult,
    TimelineResult,
    CompareResult,
]


def main() -> None:
    schema = TypeAdapter(OutputModel).json_schema()
    schema["title"] = "MolInspect output schemas"
    schema["description"] = "JSON schemas for MolInspect public API result models."

    output_path = Path(__file__).resolve().parents[1] / "schemas" / "output_schemas.json"
    output_path.write_text(json.dumps(schema, indent=2, sort_keys=True) + "\n")


if __name__ == "__main__":
    main()
