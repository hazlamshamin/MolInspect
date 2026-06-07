from __future__ import annotations

from molinspect import load
from molinspect.backends.registry import (
    OPTIONAL_BACKEND_ANNOTATIONS,
    POCKET_BACKEND_ANNOTATIONS,
    validate_backend_registry,
)
from molinspect.definitions import (
    CONTACT_CUTOFF_A,
    REFERENCE_CATALOG,
    STRUCTURAL_DEFINITIONS,
    definitions_by_kind,
    validate_definition_registry,
)
from molinspect.inspection.scales import context_scale_specs


def test_structural_definition_registry_is_internally_consistent():
    validate_definition_registry()

    assert "hydrogen_bond" in STRUCTURAL_DEFINITIONS
    assert "plip_docs" in REFERENCE_CATALOG
    assert len(STRUCTURAL_DEFINITIONS) == len(set(STRUCTURAL_DEFINITIONS))

    relation_definitions = definitions_by_kind("relation")
    assert relation_definitions
    assert all(definition.backend for definition in relation_definitions)
    assert all(definition.category for definition in relation_definitions)
    assert all(definition.default_confidence for definition in relation_definitions)


def test_backend_registry_maps_to_existing_structural_definitions():
    validate_backend_registry()

    definition_ids = {
        definition_id
        for backend_annotation in (*OPTIONAL_BACKEND_ANNOTATIONS, *POCKET_BACKEND_ANNOTATIONS)
        for definition_id in backend_annotation.definition_ids
    }
    assert "biological_interface" in definition_ids
    assert "pocket" in definition_ids
    assert definition_ids.issubset(STRUCTURAL_DEFINITIONS)


def test_contact_cutoff_is_owned_by_structural_definition_registry():
    contact = STRUCTURAL_DEFINITIONS["nonbonded_contact"]

    assert contact.parameters["max_distance_A"] == CONTACT_CUTOFF_A == 4.0
    assert contact.source == "literature_informed_heuristic"


def test_context_scales_are_owned_by_structural_definition_registry():
    context_scale = STRUCTURAL_DEFINITIONS["context_scale"]
    specs = {spec.id: spec for spec in context_scale_specs()}

    assert context_scale.source == "molinspect_policy"
    assert set(specs) == set(context_scale.parameters)
    assert specs["ligand_binding_site"].radius_A == 4.0
    assert specs["ligand_binding_site"].budget == "medium"
    assert specs["ligand_binding_site"].focus == ("ligand_contact_shell",)


def test_relation_outputs_carry_definition_attribution(tiny_pdb):
    session = load(structure=tiny_pdb)

    result = session.context("chain A and resid 2", radius=3.0, focus="contacts")
    relation = result.relations[0]

    assert relation.definition_id == relation.type
    assert relation.definition_source
    assert relation.backend
    assert relation.cutoff_A is not None
    assert relation.reference_keys


def test_structural_profiles_carry_definition_attribution(tiny_pdb):
    session = load(structure=tiny_pdb)

    result = session.locate("chain A and resid 1")
    profile = result.location.structural_profile

    assert profile is not None
    assert "local_packing_density" in profile.definition_ids
    assert profile.exposure_source is None or profile.exposure_source == "freesasa"
    assert result.location.objects == ["residue:A:1:ALA"]
