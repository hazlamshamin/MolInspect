"""Hidden structural-biology backend integrations.

Backends are implementation details behind MolInspect's small public API. Each
adapter maps external tool output onto the structural definitions registry
before domain modules consume it.
"""

from .registry import available_backend_annotation_ids

__all__ = ["available_backend_annotation_ids"]
