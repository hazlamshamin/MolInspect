from __future__ import annotations

import pytest

from molinspect import load, timeline_metrics
from molinspect.errors import MetricError


def test_distance_timeline_reports_summary_and_contact_event(tiny_pdb, tiny_multiframe_pdb):
    session = load(topology=tiny_pdb, trajectory=tiny_multiframe_pdb)

    result = session.timeline(
        metric="distance",
        selection1="chain A and resid 2",
        selection2="resname ATP",
    )

    assert result.frames_analyzed == 2
    assert result.summary["min_A"] == 2.5
    assert result.summary["max_A"] == 7.5
    assert result.summary["contact_occupancy"] == 0.5
    assert result.summary_text == (
        "distance over 2 frames: min 2.5 A, median 5 A, max 7.5 A; "
        "contact occupancy 0.5; dominant relation steric_clash."
    )
    assert {"type": "contact_breaks", "frame": 1, "distance_A": 7.5, "cutoff_A": 4.0} in result.events
    assert any(event["type"] == "ligand_moved_away" for event in result.events)
    assert result.sampled_values == [{"frame": 0, "value_A": 2.5}, {"frame": 1, "value_A": 7.5}]


def test_contact_timeline_reports_boolean_samples(tiny_pdb, tiny_multiframe_pdb):
    session = load(topology=tiny_pdb, trajectory=tiny_multiframe_pdb)

    result = session.timeline(
        metric="contact",
        selection1="chain A and resid 2",
        selection2="resname ATP",
    )

    assert result.sampled_values == [{"frame": 0, "contact": True}, {"frame": 1, "contact": False}]
    assert result.summary["contact_cutoff_A"] == 4.0
    assert result.summary["relation_type_counts"] == {"near": 1, "steric_clash": 1}
    assert result.summary["relation_occupancy"] == {"near": 0.5, "steric_clash": 0.5}


def test_relation_timeline_reports_atom_pair_evidence(tiny_pdb, tiny_multiframe_pdb):
    session = load(topology=tiny_pdb, trajectory=tiny_multiframe_pdb)

    result = session.timeline(
        metric="relation",
        selection1="chain A and resid 2",
        selection2="resname ATP",
    )

    assert result.summary["dominant_relation_type"] in {"steric_clash", "near"}
    assert result.summary["relation_type_counts"] == {"near": 1, "steric_clash": 1}
    assert result.summary["representative_frames"]["closest"]["frame"] == 0
    assert result.summary["representative_frames"]["farthest"]["frame"] == 1
    assert result.sampled_values == [
        {
            "frame": 0,
            "relation_type": "steric_clash",
            "distance_A": 2.5,
            "source_atom": "A:2:GLY:C",
            "target_atom": "A:101:ATP:P",
        },
        {
            "frame": 1,
            "relation_type": "near",
            "distance_A": 7.5,
            "source_atom": "A:2:GLY:C",
            "target_atom": "A:101:ATP:P",
        },
    ]


def test_compare_reports_lost_nearby_ligand(tiny_pdb, tiny_multiframe_pdb):
    session = load(topology=tiny_pdb, trajectory=tiny_multiframe_pdb)

    result = session.compare("chain A and resid 2", frame_a=0, frame_b=-1, radius=3.0)

    assert result.frame_a == 0
    assert result.frame_b == 1
    assert result.summary == "0 nearby objects gained, 1 lost, 1 shared within 3 A; selection RMSD 0.0 A."
    assert {"type": "nearby_object_loss", "object": "ligand:A:101:ATP", "frame": 1} in result.main_changes
    assert any(change["type"] == "selection_rmsd" for change in result.main_changes)


def test_timeline_metrics_are_public_and_concrete(tiny_pdb):
    session = load(structure=tiny_pdb)

    result = timeline_metrics()
    session_result = session.timeline_metrics()

    assert session_result == result
    assert result.count == 13
    by_metric = {metric.metric: metric for metric in result.metrics}
    assert by_metric["interaction_persistence"].required_arguments == ["selection1", "selection2"]
    assert by_metric["states"].required_arguments == ["selection"]
    assert "contact_occupancy" in by_metric["interaction_persistence"].summary_keys
    assert "state_transition" in by_metric["states"].event_types
    assert "session.timeline" in by_metric["ligand_stability"].example
    assert "not natural-language" in result.limitations[0]


def test_timeline_metric_error_explains_supported_notation(tiny_pdb, tiny_multiframe_pdb):
    session = load(topology=tiny_pdb, trajectory=tiny_multiframe_pdb)

    with pytest.raises(MetricError, match="distance.*contact.*relation.*rmsd.*rmsf.*mobility"):
        session.timeline(metric="torsion", selection1="chain A", selection2="ligand")


def test_timeline_frame_range_error_explains_notation(tiny_pdb, tiny_multiframe_pdb):
    session = load(topology=tiny_pdb, trajectory=tiny_multiframe_pdb)

    with pytest.raises(MetricError, match="inclusive .*start, end"):
        session.timeline(
            metric="distance",
            selection1="chain A",
            selection2="ligand",
            frames=[0, 1],
        )
