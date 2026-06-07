# Usage

MolInspect exposes one session object per loaded structure or trajectory. Start
with `load()`, then ask for objects, locations, local context, timelines, or
frame comparisons.

## Load A Structure

```python
from molinspect import load

session = load(structure="examples/data/4hhb.pdb")
summary = session.summary()
```

For trajectories:

```python
session = load(topology="protein.gro", trajectory="traj.xtc")
```

## Discover Objects

```python
session.objects(type="chain")
session.objects(type="residue", limit=10)
session.objects(type=["ligand", "ion"])
session.objects(type="ligand_contact_shell", contains="HEM")
session.objects(type="interchain_contact_interface")
```

Object IDs returned by `objects()` can be selected again:

```python
session.resolve_selection("ligand:A:142:HEM")
session.context("ligand_contact_shell:A:142:HEM")
```

## Locate One Selection

Use `locate()` when you want a compact identity and placement card.

```python
located = session.locate("chain A and resid 87")

print(located.resolved_objects)
print(located.location.chain)
print(located.location.secondary_structure)
print(located.location.exposure_status)
print(located.location.ligand_contact_shell_ids)
print(located.location.nearest_protein_residues[:3])
```

## Retrieve Local Context

Use `context()` when you need nearby objects and typed relation evidence.

```python
result = session.context("chain A and resname HEM", scale="ligand_binding_site")

print(result.summary)
print(result.objects[:5])
print(result.relations[:5])
print(result.limitations)
```

Useful context scales:

```text
chemical_contacts
residue_environment
ligand_binding_site
metal_coordination
protein_interface
broad_environment
```

Use `context_scales()` to inspect their concrete radius, budget, and focus
settings.

## Focus Context On A Relation Family

```python
session.context("chain A and resid 87", radius=4.0, focus="contacts")
session.context("chain A and resname HEM and name FE", radius=3.0,
                focus="metal_coordination")
session.context("chain A and resid 11", radius=5.0, focus="salt_bridges")
session.context("chain A", scale="protein_interface")
```

Focus values are explicit controls, not natural-language questions. See
`docs/API_NOTATION.md` for the full list.

## Inspect Time

Timelines work on static structures too, but are most useful with trajectories
or multi-model files.

Use `timeline_metrics()` to inspect exact metric IDs, required arguments,
example calls, summary keys, sampled-value keys, event types, and limitations.

```python
from molinspect import timeline_metrics

timeline_metrics()
session.timeline_metrics()

session.timeline(metric="distance", selection1="chain A and resid 57",
                 selection2="resname ATP")
session.timeline(metric="relation", selection1="chain A and resid 57",
                 selection2="resname ATP")
session.timeline(metric="hydrogen_bonds", selection1="chain A and resid 57",
                 selection2="resname ATP")
session.timeline(metric="interaction_persistence", selection1="chain A and resid 57",
                 selection2="resname ATP")
session.timeline(metric="ligand_stability", selection1="resname ATP",
                 selection2="around 4 of resname ATP and protein")
session.timeline(metric="centroid_distance", selection1="resid 1-29",
                 selection2="resid 122-159")
session.timeline(metric="states", selection="chain A")
session.timeline(metric="rmsd", selection="chain A")
session.timeline(metric="mobility", selection="chain A")
session.timeline(metric="displacement", selection="chain A")
session.timeline(metric="selection_spread", selection="ligand_contact_shell:A:142:HEM")
```

Compare two frames:

```python
session.compare("around 8 of resname ATP", frame_a="first", frame_b="last")
```

`metric="hydrogen_bonds"` uses MDAnalysis HydrogenBondAnalysis and requires
explicit hydrogens or donor-hydrogen topology that MDAnalysis can resolve. When
that evidence is unavailable, the result reports the limitation instead of
claiming hydrogen-bond occupancy.

`metric="interaction_persistence"` combines distance/contact, relation, and
explicit-H-bond persistence evidence. `metric="ligand_stability"` adds
site-aligned ligand RMSD and a structural stability class. `metric="states"`
returns simple representative conformational states from aligned RMSD and
selection-spread features; it is not a kinetic MSM.

Every `TimelineResult` includes:

- `summary_text`: a short interpretation sentence derived from exact fields;
- `summary`: exact metric values and representative frames;
- `events`: detected changes such as contact breaks or state transitions;
- `sampled_values`: compact frame-level evidence;
- `limitations`: backend and interpretation caveats.

## Run Included Examples

```bash
uv run python examples/real_pdb_usage.py
uv run python scripts/render_example_pdbs.py
uv run python scripts/render_example_pdbs_pymol.py
uv run python scripts/check_render_artifacts.py
```
