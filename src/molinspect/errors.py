"""MolInspect exception types."""

from __future__ import annotations


class MolInspectError(Exception):
    """Base class for package-specific errors."""


class LoadError(MolInspectError):
    """Raised when a structure or trajectory cannot be loaded."""


class FrameError(MolInspectError):
    """Raised when a frame reference is invalid for the loaded target."""


class SelectionResolutionError(MolInspectError):
    """Raised when a public selection cannot be resolved by the backend."""


class ObjectQueryError(MolInspectError):
    """Raised when object listing parameters are invalid."""


class MetricError(MolInspectError):
    """Raised when a metric cannot be computed from the requested inputs."""
