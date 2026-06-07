"""Render real PDB examples with PyMOL for visual QA.

This uses the external `pymol` executable when it is available. PyMOL remains an
optional validation/rendering tool, not a required MolInspect runtime dependency.
"""

from __future__ import annotations

import shutil
import subprocess
import tempfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
DATA_DIR = ROOT / "examples" / "data"
RENDER_DIR = ROOT / "examples" / "renders"
PYMOL = shutil.which("pymol")


def main() -> None:
    if PYMOL is None:
        raise SystemExit("PyMOL executable not found on PATH.")

    RENDER_DIR.mkdir(parents=True, exist_ok=True)
    for path in sorted(DATA_DIR.glob("*.pdb")):
        output = RENDER_DIR / f"{path.stem}_pymol.png"
        _render_with_pymol(path, output)
        print(f"wrote {output}")


def _render_with_pymol(path: Path, output: Path) -> None:
    pml = _pymol_script(path, output)
    with tempfile.TemporaryDirectory() as tmpdir:
        script_path = Path(tmpdir) / "render.pml"
        script_path.write_text(pml)
        subprocess.run([PYMOL, "-cq", str(script_path)], check=True)


def _pymol_script(path: Path, output: Path) -> str:
    pdb_path = path.resolve().as_posix()
    output_path = output.resolve().as_posix()
    return f"""
reinitialize
load "{pdb_path}", mol
remove solvent
hide everything, all
show cartoon, polymer.protein
color slate, polymer.protein
color firebrick, ss h and polymer.protein
color goldenrod, ss s and polymer.protein
color forest, ss l and polymer.protein
show sticks, organic
show spheres, inorganic
color gray70, elem C and organic
color red, elem O and organic
color blue, elem N and organic
color orange, elem P and organic
color firebrick, elem FE
set stick_radius, 0.16
set sphere_scale, 0.32
set cartoon_fancy_helices, on
set cartoon_smooth_loops, on
set depth_cue, off
set ray_opaque_background, on
set antialias, 2
set orthoscopic, on
bg_color white
orient all
zoom all, 8
clip slab, 200
png {output_path}, 1600, 1200, dpi=200, ray=1
quit
"""


if __name__ == "__main__":
    main()
