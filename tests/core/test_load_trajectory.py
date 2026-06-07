from __future__ import annotations

from molinspect import load
from molinspect.errors import FrameError


def test_load_multiframe_pdb_as_trajectory(tiny_pdb, tiny_multiframe_pdb):
    session = load(topology=tiny_pdb, trajectory=tiny_multiframe_pdb)

    assert session.n_frames == 2
    assert session.summary().mode == "trajectory"
    assert session.normalize_frame("first") == 0
    assert session.normalize_frame("last") == 1
    assert session.normalize_frame(-1) == 1


def test_invalid_frame_raises(tiny_pdb):
    session = load(structure=tiny_pdb)

    try:
        session.normalize_frame(1)
    except FrameError as exc:
        assert "outside the valid range" in str(exc)
    else:
        raise AssertionError("Expected FrameError")


def test_invalid_frame_error_explains_notation(tiny_pdb):
    session = load(structure=tiny_pdb)

    try:
        session.normalize_frame("middle")
    except FrameError as exc:
        assert "zero-based integer frame index" in str(exc)
        assert "'first'" in str(exc)
        assert "'last'" in str(exc)
    else:
        raise AssertionError("Expected FrameError")
