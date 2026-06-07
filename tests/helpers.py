from __future__ import annotations

from pathlib import Path


def repo_root() -> Path:
    return Path(__file__).resolve().parents[1]


def example_data_dir() -> Path:
    return repo_root() / "examples" / "data"


def pdb_atom_line(
    record: str,
    serial: int,
    name: str,
    resname: str,
    chain: str,
    resid: int,
    x: float,
    y: float,
    z: float,
    element: str,
) -> str:
    return (
        f"{record:<6}{serial:5d} {name:^4} {resname:>3} {chain:1}{resid:4d}    "
        f"{x:8.3f}{y:8.3f}{z:8.3f}{1.00:6.2f}{20.00:6.2f}          {element:>2}\n"
    )


def tiny_pdb_text(offset: float = 0.0) -> str:
    atoms = [
        ("ATOM", 1, "N", "ALA", "A", 1, 0.0, 0.0, 0.0, "N"),
        ("ATOM", 2, "CA", "ALA", "A", 1, 1.5, 0.0, 0.0, "C"),
        ("ATOM", 3, "C", "ALA", "A", 1, 2.5, 0.0, 0.0, "C"),
        ("ATOM", 4, "N", "GLY", "A", 2, 4.0, 0.0, 0.0, "N"),
        ("ATOM", 5, "CA", "GLY", "A", 2, 5.5, 0.0, 0.0, "C"),
        ("ATOM", 6, "C", "GLY", "A", 2, 6.5, 0.0, 0.0, "C"),
        ("HETATM", 7, "P", "ATP", "A", 101, 9.0 + offset, 0.0, 0.0, "P"),
        ("HETATM", 8, "O1", "ATP", "A", 101, 10.5 + offset, 0.0, 0.0, "O"),
        ("HETATM", 9, "O2", "ATP", "A", 101, 11.5 + offset, 0.0, 0.0, "O"),
        ("HETATM", 10, "O", "HOH", "A", 201, 20.0, 0.0, 0.0, "O"),
        ("HETATM", 11, "NA", "NA", "A", 301, 22.0, 0.0, 0.0, "NA"),
    ]
    return "".join(pdb_atom_line(*atom) for atom in atoms) + "END\n"
