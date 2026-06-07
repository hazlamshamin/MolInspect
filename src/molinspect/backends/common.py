"""Shared helpers for optional backend adapters."""

from __future__ import annotations

import importlib.util
import sys
from pathlib import Path

STATIC_STRUCTURE_SUFFIXES = frozenset({".pdb", ".ent", ".cif", ".mmcif"})


def module_is_available(module_name: str) -> bool:
    """Return whether a Python backend module is importable."""

    try:
        return importlib.util.find_spec(module_name) is not None
    except ValueError:
        return module_name in sys.modules


def module_can_import(module_name: str) -> bool:
    """Return whether a backend module imports without raising."""

    try:
        __import__(module_name)
    except Exception:
        return False
    return True


def static_structure_source(source_files: tuple[Path, ...]) -> Path | None:
    """Return the first static structure source usable by backend tools."""

    for path in source_files:
        if path.suffix.lower() in STATIC_STRUCTURE_SUFFIXES:
            return path
    return None


def has_static_structure(source_files: tuple[Path, ...]) -> bool:
    """Return whether the source set contains a static PDB/mmCIF-like file."""

    return static_structure_source(source_files) is not None


def source_cache_key(source: Path) -> tuple[str, int, int]:
    """Return a cache key that changes when a source file changes."""

    stat = source.stat()
    return (str(source.resolve()), stat.st_mtime_ns, stat.st_size)


def short_process_error(output: str) -> str:
    """Return a compact stderr/stdout tail for backend limitations."""

    lines = [line.strip() for line in output.splitlines() if line.strip()]
    if not lines:
        return "process exited with a non-zero status"
    return " | ".join(lines[-3:])[:600]
