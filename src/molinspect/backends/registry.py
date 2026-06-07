"""Central registry for optional backend availability.

This module keeps load-time annotation reporting aligned with the backend
adapters that actually provide those annotations. It is the single place where
external tools are translated into MolInspect annotation IDs.
"""

from __future__ import annotations

import shutil
from collections.abc import Callable
from dataclasses import dataclass
from pathlib import Path

from .common import has_static_structure, module_is_available
from .interfaces import pisa_is_available
from .interactions import interaction_backend_is_available, plip_is_available
from .pockets import fpocket_is_available, p2rank_is_available
from ..definitions import STRUCTURAL_DEFINITIONS


@dataclass(frozen=True, slots=True)
class BackendAnnotation:
    """One load-summary annotation produced by an optional backend."""

    annotation_id: str
    definition_ids: tuple[str, ...]
    is_available: Callable[[], bool]


BASE_ANNOTATION_IDS: tuple[str, ...] = (
    "basic_topology",
    "local_packing_density",
    "interchain_contact_interface",
    "ligand_contact_shell",
)

OPTIONAL_BACKEND_ANNOTATIONS: tuple[BackendAnnotation, ...] = (
    BackendAnnotation(
        annotation_id="dssp_secondary_structure",
        definition_ids=("secondary_structure", "loop"),
        is_available=lambda: shutil.which("mkdssp") is not None,
    ),
    BackendAnnotation(
        annotation_id="freesasa_exposure",
        definition_ids=("freesasa_exposure",),
        is_available=lambda: module_is_available("freesasa"),
    ),
    BackendAnnotation(
        annotation_id="pdbe_arpeggio_interactions",
        definition_ids=(
            "metal_coordination",
            "salt_bridge",
            "hydrogen_bond",
            "polar_contact_candidate",
            "hydrophobic_contact",
            "pi_stacking",
            "steric_clash",
            "nonbonded_contact",
        ),
        is_available=interaction_backend_is_available,
    ),
    BackendAnnotation(
        annotation_id="plip_interactions",
        definition_ids=(
            "metal_coordination",
            "salt_bridge",
            "hydrogen_bond",
            "hydrophobic_contact",
            "pi_stacking",
        ),
        is_available=plip_is_available,
    ),
    BackendAnnotation(
        annotation_id="pisa_biological_interfaces",
        definition_ids=("biological_interface",),
        is_available=pisa_is_available,
    ),
)

POCKET_BACKEND_ANNOTATIONS: tuple[BackendAnnotation, ...] = (
    BackendAnnotation(
        annotation_id="p2rank_pockets",
        definition_ids=("pocket",),
        is_available=p2rank_is_available,
    ),
    BackendAnnotation(
        annotation_id="fpocket_pockets",
        definition_ids=("pocket",),
        is_available=fpocket_is_available,
    ),
)


def available_backend_annotation_ids(source_files: tuple[Path, ...]) -> list[str]:
    """Return load-summary annotation IDs from registered backend availability."""

    annotations = list(BASE_ANNOTATION_IDS)
    if not has_static_structure(source_files):
        return annotations

    for backend_annotation in OPTIONAL_BACKEND_ANNOTATIONS:
        if backend_annotation.is_available():
            annotations.append(backend_annotation.annotation_id)

    for pocket_backend_annotation in POCKET_BACKEND_ANNOTATIONS:
        if pocket_backend_annotation.is_available():
            annotations.append(pocket_backend_annotation.annotation_id)
            break

    return annotations


def validate_backend_registry() -> None:
    """Validate backend annotation mappings against structural definitions."""

    for backend_annotation in (*OPTIONAL_BACKEND_ANNOTATIONS, *POCKET_BACKEND_ANNOTATIONS):
        for definition_id in backend_annotation.definition_ids:
            if definition_id not in STRUCTURAL_DEFINITIONS:
                raise ValueError(
                    f"{backend_annotation.annotation_id!r} maps unknown definition {definition_id!r}"
                )


validate_backend_registry()
