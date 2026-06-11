# MolInspect: Python Toolkit for Protein Structure and Trajectory Inspection

MolInspect is a Python toolkit for protein structure analysis and molecular
trajectory inspection. It turns PDB, mmCIF, and trajectory-like inputs into
compact typed evidence for structural biology, protein engineering,
ligand-binding-site analysis, protein-protein interface analysis, and LLM-agent
workflows.

A static structure is treated as a one-frame trajectory, so the same API works
for both:

```text
T = 1  -> static structure
T = N  -> trajectory, ensemble, or molecular movie
```

## Capabilities

MolInspect is useful when a workflow needs compact structural evidence instead
of raw coordinate dumps:

- Load PDB, mmCIF, and trajectory-like molecular data through MDAnalysis.
- Inspect residues, ligands, ions, waters, chains, secondary-structure elements,
  ligand contact shells, pockets, and interfaces.
- Retrieve residue environments, ligand-binding-site context, metal-coordination
  context, and protein-interface context with explicit selection notation.
- Report typed interaction evidence for hydrogen bonds, salt bridges, metal
  coordination, hydrophobic contacts, pi stacking, water-bridge candidates,
  steric clashes, nonbonded contacts, and nearby residues.
- Summarize trajectory and ensemble signals: distance, contact occupancy,
  relation occupancy, hydrogen-bond occupancy, ligand stability, RMSD,
  RMSF/mobility, displacement, conformational states, centroid distance, and
  selection spread.
- Return Pydantic models with exact measurements, object IDs, source atoms,
  backend attribution, cutoffs, references, and limitations.

## Who It Is For

- Structural biologists inspecting local protein structure context.
- Protein engineers checking mutation neighborhoods, ligand-facing residues, or
  interface-facing residues.
- Computational biologists and bioinformatics workflows that need structured
  PDB/mmCIF evidence.
- Molecular dynamics users who need compact trajectory summaries for RMSD, RMSF,
  contact persistence, ligand stability, and interaction changes.
- LLM agents and analysis pipelines that need high-signal structural context
  without reading entire coordinate files.

## Backends And Evidence

MolInspect uses a small public API over internal structural-biology backends:

- MDAnalysis for molecular loading, selection, distance, frame, RMSD, RMSF, and
  hydrogen-bond trajectory support.
- FreeSASA for exact solvent-accessible surface area and exposure.
- DSSP through `mkdssp` for secondary-structure and loop objects.
- PLIP and PDBe Arpeggio for protein-ligand and non-covalent interaction
  evidence.
- PISA for biological-interface objects.
- P2Rank for ligand-binding-site pocket prediction; fpocket for geometric
  pocket detection.

It is not a molecular simulator, docking engine, protein-design tool, or full
visualization GUI. It does not estimate binding free energy or kinetic state
models.

## Install

```bash
uv sync
```

For local development, use:

```bash
uv sync --extra dev
```

For the full backend stack, put these commands on `PATH`: `mkdssp` for DSSP
secondary structure, `pisa` for biological interfaces, and `prank`/`p2rank` or
`fpocket` for pocket objects. mmCIF loading uses the Gemmi extra:

```bash
uv sync --extra static
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
session.summary()
session.resolve_selection()
session.objects()
session.locate()
session.context()
session.context_scales()
session.timeline_metrics()
session.timeline()
session.compare()
```

Typical calls:

```python
session.objects(type="ligand")
session.resolve_selection("ligand:A:142:HEM")
session.locate("chain A and resid 87")
session.context("chain A and resid 87", scale="residue_environment")
session.context("chain A and resname HEM", scale="ligand_binding_site")
session.context("protein and chain A", scale="protein_interface")
session.context_scales()
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
- `llms.txt` lists the canonical docs for LLM and AI-search ingestion.

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
uv run mypy src tests
uv build
```
