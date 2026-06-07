from __future__ import annotations

from pathlib import Path

import pytest

from molinspect import load
from molinspect.errors import MetricError


def _pdb_atom_line(
    serial: int,
    name: str,
    x: float,
    y: float,
    z: float,
) -> str:
    return (
        f"{'ATOM':<6}{serial:5d} {name:^4} {'ALA':>3} {'A':1}{1:4d}    "
        f"{x:8.3f}{y:8.3f}{z:8.3f}{1.00:6.2f}{20.00:6.2f}          {name[0]:>2}\n"
    )


def _frame_text(ca_y: float, c_y: float) -> str:
    atoms = [
        _pdb_atom_line(1, "N", 0.0, 0.0, 0.0),
        _pdb_atom_line(2, "CA", 1.4, ca_y, 0.0),
        _pdb_atom_line(3, "C", 2.8, c_y, 0.0),
        _pdb_atom_line(4, "O", 3.7, c_y, 0.0),
    ]
    return "".join(atoms)


def _deforming_trajectory(tmp_path: Path) -> tuple[Path, Path]:
    topology = tmp_path / "deforming_topology.pdb"
    trajectory = tmp_path / "deforming_trajectory.pdb"
    frame_0 = _frame_text(ca_y=0.0, c_y=0.0)
    frame_1 = _frame_text(ca_y=0.2, c_y=0.9)
    frame_2 = _frame_text(ca_y=0.8, c_y=4.0)
    topology.write_text(frame_0 + "END\n")
    trajectory.write_text(
        f"MODEL        1\n{frame_0}ENDMDL\n"
        f"MODEL        2\n{frame_1}ENDMDL\n"
        f"MODEL        3\n{frame_2}ENDMDL\nEND\n"
    )
    return topology, trajectory


def _hbond_atom_line(
    serial: int,
    name: str,
    resname: str,
    resid: int,
    y: float,
    element: str,
) -> str:
    return (
        f"{'ATOM':<6}{serial:5d} {name:^4} {resname:>3} {'A':1}{resid:4d}    "
        f"{0.0:8.3f}{y:8.3f}{0.0:8.3f}{1.00:6.2f}{20.00:6.2f}          {element:>2}\n"
    )


def _hbond_frame_text(acceptor_y: float) -> str:
    atoms = [
        _hbond_atom_line(1, "N", "ALA", 1, 0.0, "N"),
        _hbond_atom_line(2, "H", "ALA", 1, 1.0, "H"),
        _hbond_atom_line(3, "O", "GLY", 2, acceptor_y, "O"),
    ]
    return "".join(atoms)


def _hydrogen_bond_trajectory(tmp_path: Path) -> tuple[Path, Path]:
    topology = tmp_path / "hbond_topology.pdb"
    trajectory = tmp_path / "hbond_trajectory.pdb"
    frame_0 = _hbond_frame_text(acceptor_y=2.7)
    frame_1 = _hbond_frame_text(acceptor_y=5.0)
    topology.write_text(frame_0 + "END\n")
    trajectory.write_text(
        f"MODEL        1\n{frame_0}ENDMDL\n"
        f"MODEL        2\n{frame_1}ENDMDL\nEND\n"
    )
    return topology, trajectory


def _large_state_trajectory(tmp_path: Path) -> tuple[Path, Path]:
    topology = tmp_path / "state_topology.pdb"
    trajectory = tmp_path / "state_trajectory.pdb"
    frame_0 = _frame_text(ca_y=0.0, c_y=0.0)
    frame_1 = _frame_text(ca_y=0.2, c_y=0.9)
    frame_2 = _frame_text(ca_y=2.5, c_y=8.0)
    topology.write_text(frame_0 + "END\n")
    trajectory.write_text(
        f"MODEL        1\n{frame_0}ENDMDL\n"
        f"MODEL        2\n{frame_1}ENDMDL\n"
        f"MODEL        3\n{frame_2}ENDMDL\nEND\n"
    )
    return topology, trajectory


def test_rmsd_timeline_reports_alignment_and_representative_frames(tmp_path):
    topology, trajectory = _deforming_trajectory(tmp_path)
    session = load(topology=topology, trajectory=trajectory)

    result = session.timeline(metric="rmsd", selection="chain A and resid 1")

    assert result.frames_analyzed == 3
    assert result.summary["reference_frame"] == 0
    assert result.summary["alignment"] == "kabsch_over_selection"
    assert result.summary["min_rmsd_A"] == 0.0
    assert result.summary["max_rmsd_A"] > 0.1
    assert result.summary["representative_frames"]["highest"]["frame"] == 2
    assert result.sampled_values[0] == {"frame": 0, "rmsd_A": 0.0}
    assert result.events[0]["type"] == "rmsd_peak"


def test_mobility_timeline_reports_residue_rmsf(tmp_path):
    topology, trajectory = _deforming_trajectory(tmp_path)
    session = load(topology=topology, trajectory=trajectory)

    result = session.timeline(metric="mobility", selection="chain A and resid 1")

    sample = result.sampled_values[0]
    assert result.summary["max_rmsf_A"] > 0.1
    assert result.summary["most_mobile_object"] == "residue:A:1:ALA"
    assert sample["object"] == "residue:A:1:ALA"
    assert sample["mean_rmsf_A"] > 0.0
    assert sample["max_atom_rmsf_A"] >= sample["mean_rmsf_A"]
    assert sample["selected_atom_count"] == 4
    assert result.events[0]["type"] == "mobility_peak"


def test_hydrogen_bond_timeline_reports_backend_occupancy(tmp_path):
    topology, trajectory = _hydrogen_bond_trajectory(tmp_path)
    session = load(topology=topology, trajectory=trajectory)

    result = session.timeline(
        metric="hydrogen_bonds",
        selection1="chain A and resid 1",
        selection2="chain A and resid 2",
    )

    assert result.summary["backend"] == "MDAnalysis HydrogenBondAnalysis"
    assert result.summary["hydrogen_bond_occupancy"] == 0.5
    assert result.summary["total_hydrogen_bond_observations"] == 1
    assert result.summary["top_hydrogen_bond_pairs"][0]["occupancy"] == 0.5
    assert result.summary["representative_frames"]["first_present"]["frame"] == 0
    assert result.events == [
        {
            "type": "hydrogen_bond_breaks",
            "frame": 1,
            "hydrogen_bond_count": 0,
            "source": "residue:A:1:ALA",
            "target": "residue:A:2:GLY",
        }
    ]
    assert result.sampled_values == [
        {
            "frame": 0,
            "hydrogen_bond_count": 1,
            "has_hydrogen_bond": True,
            "donor_atom": "A:1:ALA:N",
            "hydrogen_atom": "A:1:ALA:H",
            "acceptor_atom": "A:2:GLY:O",
            "distance_A": 2.7,
            "angle_deg": 180.0,
        },
        {"frame": 1, "hydrogen_bond_count": 0, "has_hydrogen_bond": False},
    ]


def test_relation_timeline_prefers_backend_hydrogen_bond_when_available(tmp_path):
    topology, trajectory = _hydrogen_bond_trajectory(tmp_path)
    session = load(topology=topology, trajectory=trajectory)

    result = session.timeline(
        metric="relation",
        selection1="chain A and resid 1",
        selection2="chain A and resid 2",
    )

    assert result.sampled_values[0]["relation_type"] == "hydrogen_bond"
    assert result.summary["relation_occupancy"]["hydrogen_bond"] == 0.5
    assert "MDAnalysis HydrogenBondAnalysis" in result.limitations[2]


def test_hydrogen_bond_timeline_reports_missing_hydrogen_limitation(tiny_pdb, tiny_multiframe_pdb):
    session = load(topology=tiny_pdb, trajectory=tiny_multiframe_pdb)

    result = session.timeline(
        metric="hydrogen_bonds",
        selection1="chain A and resid 1",
        selection2="chain A and resid 2",
    )

    assert result.summary["hydrogen_bond_occupancy"] == 0.0
    assert result.sampled_values == [
        {"frame": 0, "hydrogen_bond_count": 0, "has_hydrogen_bond": False},
        {"frame": 1, "hydrogen_bond_count": 0, "has_hydrogen_bond": False},
    ]
    assert any("explicit hydrogens" in limitation for limitation in result.limitations)


def test_interaction_persistence_combines_contact_relation_and_hbond_evidence(tmp_path):
    topology, trajectory = _hydrogen_bond_trajectory(tmp_path)
    session = load(topology=topology, trajectory=trajectory)

    result = session.timeline(
        metric="interaction_persistence",
        selection1="chain A and resid 1",
        selection2="chain A and resid 2",
    )

    assert result.summary["contact_occupancy"] == 0.5
    assert result.summary["hydrogen_bond_occupancy"] == 0.5
    assert result.summary["relation_occupancy"] == {"hydrogen_bond": 0.5, "near": 0.5}
    assert result.summary_text == (
        "interaction_persistence over 2 frames: contact occupancy 0.5; "
        "H-bond occupancy 0.5; dominant relation hydrogen_bond."
    )
    assert result.summary["top_hydrogen_bond_pairs"][0]["occupancy"] == 0.5
    assert any(event["type"] == "hydrogen_bond_breaks" for event in result.events)
    assert result.sampled_values[0]["contact"] is True
    assert result.sampled_values[0]["has_hydrogen_bond"] is True


def test_ligand_stability_reports_site_aligned_repositioning(tiny_pdb, tiny_multiframe_pdb):
    session = load(topology=tiny_pdb, trajectory=tiny_multiframe_pdb)

    result = session.timeline(
        metric="ligand_stability",
        selection1="resname ATP",
        selection2="chain A and resid 2",
    )

    assert result.summary["alignment"] == "binding_site_kabsch"
    assert result.summary["contact_occupancy"] == 0.5
    assert result.summary["max_site_aligned_ligand_rmsd_A"] == 5.0
    assert result.summary["stability_class"] == "unstable_or_repositioned"
    assert result.summary_text == (
        "ligand_stability over 2 frames: unstable_or_repositioned; "
        "max site-aligned ligand RMSD 5 A; contact occupancy 0.5."
    )
    assert result.summary["representative_frames"]["least_stable_pose"]["frame"] == 1
    assert any(event["type"] == "ligand_pose_repositioned" for event in result.events)


def test_states_timeline_reports_simple_representative_states(tmp_path):
    topology, trajectory = _large_state_trajectory(tmp_path)
    session = load(topology=topology, trajectory=trajectory)

    result = session.timeline(metric="states", selection="chain A and resid 1")

    assert result.summary["state_count"] == 2
    assert result.summary["state_method"] == "deterministic_rmsd_split"
    assert result.summary_text.startswith(
        "states over 3 frames: 2 representative state(s) by deterministic_rmsd_split;"
    )
    assert [state["state_id"] for state in result.summary["states"]] == ["state_1", "state_2"]
    assert result.events == [
        {"type": "state_transition", "frame": 2, "from_state": "state_1", "to_state": "state_2"}
    ]
    assert result.sampled_values[-1]["state_id"] == "state_2"


def test_displacement_timeline_reports_per_residue_motion(tmp_path):
    topology, trajectory = _deforming_trajectory(tmp_path)
    session = load(topology=topology, trajectory=trajectory)

    result = session.timeline(metric="displacement", selection="chain A and resid 1")

    assert result.summary["reference_frame"] == 0
    assert result.summary["max_displacement_A"] > 0.1
    assert result.summary["most_displaced_object"] == "residue:A:1:ALA"
    assert result.summary["representative_frames"]["highest"]["frame"] == 2
    assert result.sampled_values[-1]["most_displaced_object"] == "residue:A:1:ALA"
    assert result.events[0]["type"] == "displacement_peak"


def test_selection_spread_timeline_reports_least_and_most_spread_frames(tmp_path):
    topology, trajectory = _deforming_trajectory(tmp_path)
    session = load(topology=topology, trajectory=trajectory)

    result = session.timeline(metric="selection_spread", selection="chain A and resid 1")

    assert result.summary["spread_proxy"] == "selection_radius_of_gyration"
    assert result.summary["min_selection_spread_A"] < result.summary["max_selection_spread_A"]
    assert result.summary["representative_frames"]["least_spread"]["frame"] == 0
    assert result.summary["representative_frames"]["most_spread"]["frame"] == 2
    assert any(event["type"] == "selection_spread_increased" for event in result.events)


def test_temporal_metrics_respect_frame_ranges_and_compare_reports_selection_rmsd(tmp_path):
    topology, trajectory = _deforming_trajectory(tmp_path)
    session = load(topology=topology, trajectory=trajectory)

    ranged = session.timeline(
        metric="rmsd",
        selection="chain A and resid 1",
        frames=("first", "last"),
        stride=2,
    )
    compared = session.compare("chain A and resid 1", frame_a=0, frame_b=-1, radius=4.0)

    assert ranged.frames_analyzed == 2
    assert [sample["frame"] for sample in ranged.sampled_values] == [0, 2]
    rmsd_changes = [change for change in compared.main_changes if change["type"] == "selection_rmsd"]
    assert rmsd_changes and rmsd_changes[0]["rmsd_A"] > 0.1
    assert compared.evidence[-1].metric == "selection_rmsd"


def test_ligand_motion_events_and_coordinate_anti_shortcut(tiny_pdb, tiny_multiframe_pdb):
    session = load(topology=tiny_pdb, trajectory=tiny_multiframe_pdb)

    result = session.timeline(
        metric="distance",
        selection1="chain A and resid 2",
        selection2="resname ATP",
    )

    assert result.sampled_values == [{"frame": 0, "value_A": 2.5}, {"frame": 1, "value_A": 7.5}]
    assert any(event["type"] == "ligand_moved_away" and event["frame"] == 1 for event in result.events)


def test_centroid_distance_timeline_reports_region_opening(tiny_pdb, tiny_multiframe_pdb):
    session = load(topology=tiny_pdb, trajectory=tiny_multiframe_pdb)

    result = session.timeline(
        metric="centroid_distance",
        selection1="chain A and resid 2",
        selection2="resname ATP",
    )

    assert result.summary["center"] == "center_of_geometry"
    assert result.summary["min_centroid_distance_A"] == 5.0
    assert result.summary["max_centroid_distance_A"] == 10.0
    assert result.summary["representative_frames"]["closest_centroids"] == {
        "frame": 0,
        "centroid_distance_A": 5.0,
    }
    assert result.summary["representative_frames"]["farthest_centroids"] == {
        "frame": 1,
        "centroid_distance_A": 10.0,
    }
    assert result.sampled_values == [
        {"frame": 0, "centroid_distance_A": 5.0},
        {"frame": 1, "centroid_distance_A": 10.0},
    ]
    assert result.events == [
        {
            "type": "centroid_distance_increased",
            "frame": 1,
            "from_centroid_distance_A": 5.0,
            "to_centroid_distance_A": 10.0,
            "delta_A": 5.0,
        }
    ]


def test_frame_order_reversal_is_rejected(tmp_path):
    topology, trajectory = _deforming_trajectory(tmp_path)
    session = load(topology=topology, trajectory=trajectory)

    with pytest.raises(MetricError, match="frame range end must be >= start"):
        session.timeline(
            metric="rmsd",
            selection="chain A and resid 1",
            frames=("last", "first"),
        )
