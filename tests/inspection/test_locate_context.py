from __future__ import annotations

import pytest

from molinspect import context_scales, load
from molinspect.errors import MetricError


def test_locate_returns_compact_location_evidence(tiny_pdb):
    session = load(structure=tiny_pdb)

    result = session.locate("chain A and resid 2")

    assert result.resolved_objects == ["residue:A:2:GLY"]
    assert result.frame == 0
    assert result.location.chain == "A"
    assert result.location.center_of_geometry_A == [5.333, 0.0, 0.0]
    nearest_residue = result.location.nearest_protein_residues[0]
    assert nearest_residue.object_id == "residue:A:1:ALA"
    assert nearest_residue.min_distance_A == 1.5
    assert nearest_residue.annotations is not None

    nearest_ligand = result.location.near_ligands[0]
    assert nearest_ligand.object_id == "ligand:A:101:ATP"
    assert nearest_ligand.min_distance_A == 2.5
    assert result.location.structural_profile
    assert result.location.structural_profile.local_packing_source == "local_packing_density"
    assert result.location.structural_profile.exposure_source in {None, "freesasa"}


def test_context_returns_radius_relations_without_full_graph(tiny_pdb):
    session = load(structure=tiny_pdb)

    result = session.context("chain A and resid 2", radius=3.0, focus="contacts")

    object_ids = {obj.id for obj in result.objects}
    assert result.scale is None
    assert result.radius_A == 3.0
    assert result.budget == "small"
    assert result.focus == ["contacts"]
    assert object_ids == {"residue:A:2:GLY", "residue:A:1:ALA", "ligand:A:101:ATP"}
    assert len(result.relations) == 2
    assert {(rel.target, rel.type, rel.min_distance_A) for rel in result.relations} == {
        ("residue:A:1:ALA", "inferred_covalent_bond", 1.5),
        ("ligand:A:101:ATP", "steric_clash", 2.5),
    }
    assert result.summary.startswith("2 nearby objects within 3 A; 2 are contacts at <= 4 A.")
    assert "Closest target is residue:A:1:ALA at 1.5 A." in result.summary


def test_context_scale_expands_to_visible_effective_controls(tiny_pdb):
    session = load(structure=tiny_pdb)

    result = session.context("chain A and resid 2", scale="ligand_binding_site")

    assert result.scale == "ligand_binding_site"
    assert result.radius_A == 4.0
    assert result.budget == "medium"
    assert result.focus == ["ligand_contact_shell"]
    assert any(item.metric == "context_scale" and item.value == "ligand_binding_site" for item in result.evidence)
    assert any(item.metric == "budget" and item.value == "medium" for item in result.evidence)


def test_context_scales_are_public_and_concrete(tiny_pdb):
    session = load(structure=tiny_pdb)

    result = context_scales()
    session_result = session.context_scales()

    assert session_result == result
    assert result.count == 6
    scale_by_id = {scale.id: scale for scale in result.scales}
    assert scale_by_id["ligand_binding_site"].radius_A == 4.0
    assert scale_by_id["ligand_binding_site"].budget == "medium"
    assert scale_by_id["ligand_binding_site"].focus == ["ligand_contact_shell"]
    assert "retrieval presets" in result.limitations[0]


def test_context_scale_allows_explicit_overrides(tiny_pdb):
    session = load(structure=tiny_pdb)

    result = session.context(
        "chain A and resid 2",
        scale="ligand_binding_site",
        radius=6.0,
        budget="large",
        focus="contacts",
    )

    assert result.scale == "ligand_binding_site"
    assert result.radius_A == 6.0
    assert result.budget == "large"
    assert result.focus == ["contacts"]


def test_context_rejects_unknown_scale(tiny_pdb):
    session = load(structure=tiny_pdb)

    with pytest.raises(MetricError, match="Supported context scale values"):
        session.context("chain A and resid 2", scale="zoom in")


def test_context_rejects_unknown_focus(tiny_pdb):
    session = load(structure=tiny_pdb)

    with pytest.raises(MetricError, match="Unsupported context focus"):
        session.context("chain A and resid 2", focus="what is near this")
