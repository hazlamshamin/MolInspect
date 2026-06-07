"""Quick smoke usage on real PDB example files."""

from __future__ import annotations

from pathlib import Path

from molinspect import load

ROOT = Path(__file__).resolve().parents[1]
DATA_DIR = ROOT / "examples" / "data"


def main() -> None:
    for path in sorted(DATA_DIR.glob("*.pdb")):
        session = load(structure=path)
        summary = session.summary()
        chains = session.objects(type="chain", limit=20)
        ligands = session.objects(type="ligand", limit=20)
        first_chain = chains.objects[0].chain or "_"
        first_residue = session.objects(type="residue", limit=1).objects[0]
        selection = session.resolve_selection(f"chain {first_chain} and resid {first_residue.resid}")

        print(path.name)
        print(f"  atoms={summary.n_atoms} residues={summary.n_residues} frames={summary.n_frames}")
        print(f"  chains={[obj.id for obj in chains.objects]}")
        print(f"  ligands={[obj.id for obj in ligands.objects]}")
        print(f"  selection={selection.selection!r} -> {selection.resolved_objects}")

        if path.stem == "4hhb":
            heme = session.locate("chain A and resname HEM")
            heme_context = session.context("chain A and resname HEM", radius=4.0)
            print(f"  heme_location={heme.location.plain_language}")
            print(f"  heme_context={heme_context.summary}")


if __name__ == "__main__":
    main()
