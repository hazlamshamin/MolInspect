from __future__ import annotations

from pathlib import Path

import pytest

from tests.helpers import tiny_pdb_text


@pytest.fixture
def tiny_pdb(tmp_path: Path) -> Path:
    path = tmp_path / "tiny.pdb"
    path.write_text(tiny_pdb_text())
    return path


@pytest.fixture
def tiny_multiframe_pdb(tmp_path: Path) -> Path:
    path = tmp_path / "tiny_traj.pdb"
    frame_0 = tiny_pdb_text().replace("END\n", "")
    frame_1 = tiny_pdb_text(offset=5.0).replace("END\n", "")
    path.write_text(f"MODEL        1\n{frame_0}ENDMDL\nMODEL        2\n{frame_1}ENDMDL\nEND\n")
    return path
