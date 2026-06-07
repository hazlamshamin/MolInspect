"""Render simple QA screenshots for real PDB example structures.

This is not MolInspect's public `render()` API. It is a local validation helper that
uses the current API to load structures, then plots a compact C-alpha trace plus ligand
atoms so humans can inspect that the examples are sane.
"""

from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from molinspect import load

ROOT = Path(__file__).resolve().parents[1]
DATA_DIR = ROOT / "examples" / "data"
RENDER_DIR = ROOT / "examples" / "renders"


def main() -> None:
    RENDER_DIR.mkdir(parents=True, exist_ok=True)
    for path in sorted(DATA_DIR.glob("*.pdb")):
        output = RENDER_DIR / f"{path.stem}_overview.png"
        render_overview(path, output)
        print(f"wrote {output}")


def render_overview(path: Path, output: Path) -> None:
    session = load(structure=path)
    universe = session.universe
    summary = session.summary()

    protein_ca = universe.select_atoms("protein and name CA")
    if len(protein_ca) == 0:
        protein_ca = universe.select_atoms("protein")

    ligand_selection = session.resolve_selection("ligand")
    ligand_atoms = universe.select_atoms(ligand_selection.expression)

    fig = plt.figure(figsize=(8, 6), dpi=160)
    ax = fig.add_subplot(111, projection="3d")
    fig.patch.set_facecolor("white")
    ax.set_facecolor("white")

    if len(protein_ca) > 0:
        positions = protein_ca.positions
        residue_values = np.arange(len(positions))
        ax.plot(
            positions[:, 0],
            positions[:, 1],
            positions[:, 2],
            color="#2f6f9f",
            linewidth=1.3,
            alpha=0.78,
        )
        ax.scatter(
            positions[:, 0],
            positions[:, 1],
            positions[:, 2],
            c=residue_values,
            cmap="viridis",
            s=18,
            depthshade=True,
            label="protein CA",
        )

    if len(ligand_atoms) > 0:
        ligand_positions = ligand_atoms.positions
        ax.scatter(
            ligand_positions[:, 0],
            ligand_positions[:, 1],
            ligand_positions[:, 2],
            color="#c43b2f",
            s=26,
            depthshade=True,
            label="ligand atoms",
        )

    _set_equal_axes(ax)
    ax.set_title(
        f"{path.stem.upper()} | atoms={summary.n_atoms}, residues={summary.n_residues}, "
        f"chains={summary.n_chains}, ligands={session.objects(type='ligand').count}",
        pad=16,
    )
    ax.set_xlabel("X (A)")
    ax.set_ylabel("Y (A)")
    ax.set_zlabel("Z (A)")
    ax.legend(loc="upper right")
    ax.view_init(elev=20, azim=35)
    fig.tight_layout()
    fig.savefig(output)
    plt.close(fig)


def _set_equal_axes(ax) -> None:
    limits = np.array([ax.get_xlim3d(), ax.get_ylim3d(), ax.get_zlim3d()])
    centers = limits.mean(axis=1)
    radius = max((limits[:, 1] - limits[:, 0]).max() / 2.0, 1.0)
    ax.set_xlim3d([centers[0] - radius, centers[0] + radius])
    ax.set_ylim3d([centers[1] - radius, centers[1] + radius])
    ax.set_zlim3d([centers[2] - radius, centers[2] + radius])


if __name__ == "__main__":
    main()
