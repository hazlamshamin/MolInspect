"""Evidence-backed MolInspect structural vocabulary.

The registry owns public structural terms, reference keys, confidence/source
labels, and definition-derived constants. Runtime code should import semantic
cutoffs and definition metadata from this package instead of redefining them
inside backend or domain modules.
"""

from .registry import *  # noqa: F403
