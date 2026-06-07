# MolInspect

MolInspect is a Python inspection layer for molecular structures and
trajectories. It helps programs and LLM agents retrieve compact structural
evidence from PDB, mmCIF, and trajectory-like inputs without dumping raw
coordinates into context.

A static structure is treated as a one-frame trajectory, so the same API works
for both:

```text
T = 1  -> static structure
T = N  -> trajectory, ensemble, or molecular movie
```

## What It Helps With

MolInspect is useful when you need structured answers to questions such as:

- where a residue, ligand, ion, pocket, or interface sits in a structure;
- which residues form a ligand-binding environment;
- which atoms support metal coordination, salt bridges, hydrogen-bond
  candidates, hydrophobic contacts, aromatic contacts, or steric clashes;
- whether a contact, hydrogen bond, distance, RMSD, mobility signal, or selection
  spread changes over frames;
- which compact evidence slice is worth handing to an LLM, report, or downstream
  structural-biology workflow.

It is not a molecular simulator, docking engine, protein-design tool, or full
visualization GUI.

## Install

```bash
uv sync
```

## Quick Start

```python
from molinspect import load

session = load(structure="examples/data/4hhb.pdb")

print(session.summary())
print(session.objects(type=["chain", "ligand"]))

heme = session.locate("chain A and resname HEM")
print(heme.location.plain_language)

site = session.context("chain A and resname HEM", scale="ligand_binding_site")
print(site.summary)

iron = session.context(
    "chain A and resname HEM and name FE",
    radius=3.0,
    focus="metal_coordination",
)
print(iron.relations[0])
```

## Core API

MolInspect keeps the public surface intentionally small:

```python
load()
objects()
locate()
context()
timeline()
compare()
```

Typical calls:

```python
session.objects(type="ligand")
session.locate("chain A and resid 87")
session.context("chain A and resid 87", scale="residue_environment")
session.context("chain A and resname HEM", scale="ligand_binding_site")
session.context("chain A", scale="protein_interface")
session.timeline_metrics()
session.timeline(metric="centroid_distance", selection1="resid 1-29",
                 selection2="resid 122-159")
session.timeline(metric="relation", selection1="chain A and resname HEM and name FE",
                 selection2="chain A and resid 87")
session.timeline(metric="hydrogen_bonds", selection1="chain A and resid 57",
                 selection2="resname ATP")
session.timeline(metric="interaction_persistence", selection1="chain A and resid 57",
                 selection2="resname ATP")
session.timeline(metric="ligand_stability", selection1="resname ATP",
                 selection2="around 4 of resname ATP and protein")
session.timeline(metric="states", selection="chain A")
session.timeline(metric="rmsd", selection="chain A")
session.compare("chain A and resname HEM", frame_a="first", frame_b="last")
```

## Public Docs

- `docs/USAGE.md` shows common workflows.
- `docs/APPLICATIONS.md` explains where MolInspect fits in structural-biology,
  protein-engineering, MD, and LLM-agent workflows.
- `docs/API_SPEC.md` documents the public Python API.
- `docs/API_NOTATION.md` defines the explicit selection, object ID, frame, focus,
  context-scale, metric, and relation notation.

## Examples

Real PDB examples are included for local smoke checks:

```bash
uv run python examples/real_pdb_usage.py
uv run python scripts/render_example_pdbs.py
uv run python scripts/render_example_pdbs_pymol.py
uv run python scripts/check_render_artifacts.py
```

The included examples are `1CRN`, `1UBQ`, and `4HHB`.

## Development Checks

```bash
uv sync --extra dev
uv run pytest
uv run ruff check .
uv run mypy src
uv build
```
