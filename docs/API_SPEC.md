# API Specification

This document describes the public Python API.

For exact string formats used by these calls, read `docs/API_NOTATION.md`.

All string controls are explicit structural notation. They are not
natural-language commands.

Return snippets below are compact examples. Complete output field names and
types are also available in `schemas/output_schemas.json`.

## 1. `load()`

Load a static structure or trajectory.

```python
session = load(
    structure: str | Path | None = None,
    topology: str | Path | None = None,
    trajectory: str | Path | None = None,
    name: str | None = None,
)
```

### Examples

Static:

```python
session = load(structure="1abc.cif")
```

Trajectory:

```python
session = load(topology="protein.gro", trajectory="traj.xtc")
```

### Return

`load()` returns an `InspectionSession`. The session exposes the load metadata through
`session.load_result` and `session.summary()`.

```json
{
  "target_id": "sha256:...",
  "name": "protein",
  "n_atoms": 48231,
  "n_residues": 326,
  "n_chains": 2,
  "n_segments": 2,
  "n_frames": 10000,
  "mode": "trajectory",
  "source_files": ["protein.gro", "traj.xtc"],
  "available_annotations": [
    "basic_topology",
    "local_packing_density",
    "interchain_contact_interface",
    "ligand_contact_shell",
    "dssp_secondary_structure",
    "freesasa_exposure",
    "pdbe_arpeggio_interactions",
    "plip_interactions",
    "p2rank_pockets",
    "pisa_biological_interfaces"
  ],
  "limitations": []
}
```

`fpocket_pockets` appears instead of `p2rank_pockets` when only fpocket is
discoverable. `pisa_biological_interfaces` appears only when a local `pisa`
command is discoverable.

---

## 2. Session helpers

The returned `InspectionSession` also exposes small helper methods:

```python
session.summary() -> LoadResult
session.resolve_selection(selection: str, frame: int | str = 0) -> SelectionResult
session.context_scales() -> ContextScalesResult
session.timeline_metrics() -> TimelineMetricsResult
```

`summary()` returns the same load metadata as `session.load_result`.
`resolve_selection()` is useful when users or agents need to verify exactly what
a selection string resolves to before using `locate()`, `context()`, `timeline()`,
or `compare()`.

`context_scales()` and `timeline_metrics()` expose the concrete accepted scale
and metric controls, including required arguments and output keys.

---

## 3. `objects()`

List stable molecular objects with optional explicit filters.

```python
session.objects(
    type: str | list[str] | None = None,
    contains: str | None = None,
    frame: int | str | None = None,
    limit: int = 50,
)
```

`contains` is a literal, case-insensitive substring filter over object IDs and
labels. It is not natural-language input.

Supported `type` values:

```text
atom
residue
ligand
ion
water
chain
secondary_structure
loop
ligand_contact_shell
pocket
interchain_contact_interface
biological_interface
```

Plural aliases such as `ligands`, `residues`, `pockets`,
`secondary_structures`, and `pisa_interfaces` are accepted. `selection_region`
can appear in broad-selection outputs, but it is not a queryable
`objects(type=...)` value.

### Example

```python
session.objects(type=["ligand", "ion"])
```

### Object ID policy

Object IDs are stable within a loaded target and do not require callers to reason
from raw atom indices alone.

Canonical forms:

```text
atom:{chain_or_segment}:{resid}{icode}:{resname}:{atom_name}:{index0}
residue:{chain_or_segment}:{resid}{icode}:{resname}
ligand:{chain_or_segment}:{resid}{icode}:{resname}
ion:{chain_or_segment}:{resid}{icode}:{resname}
water:{chain_or_segment}:{resid}{icode}:{resname}
chain:{chain_or_segment}
selection_region:{n_residues}res:{first_atom_index0}-{last_atom_index0}:{fingerprint}
secondary_structure:{chain}:{kind}:{start_resid}-{end_resid}
loop:{chain}:{start_resid}-{end_resid}
ligand_contact_shell:{chain}:{resid}{icode}:{resname}
pocket:p2rank:{rank}
pocket:fpocket:{rank}
interchain_contact_interface:{chain_a}-{chain_b}
biological_interface:pisa:{interface_id}
```

When a chain ID is unavailable, use the segment ID. When neither is meaningful, use
`_`. The zero-based atom index is included only for atoms because atom names can repeat
within unusual residues or topology sources.

Examples:

```text
residue:A:87:HIS
ligand:A:142:HEM
atom:A:142:HEM:FE:1064
biological_interface:pisa:1
```

### Return

```json
{
  "objects": [
    {
      "id": "ligand:A:401:ATP",
      "type": "ligand",
      "name": "ATP",
      "chain": "A",
      "resid": 401,
      "atom_count": 31
    }
  ],
  "count": 1,
  "returned": 1,
  "truncated": false,
  "limitations": []
}
```

---

## 4. `locate()`

Answer “where is this object/region?”

```python
session.locate(
    selection: str,
    frame: int | str = 0,
    include_metrics: bool = True,
)
```

`frame` can be:

- integer zero-based frame index;
- negative integer frame index from the end, such as `-1`;
- `"first"`;
- `"last"`.

Unsupported frame labels raise `FrameError`.

### Selection Notation

Stable selection examples:

```text
chain A and resid 57
A:57
resname ATP
name FE
ligand
water
ion
around 4 of chain A and resid 57
```

### Return

```json
{
  "selection": "chain A and resid 57",
  "resolved_objects": ["residue:A:57:HIS"],
  "frame": 0,
  "location": {
    "center_of_geometry_A": [12.3, 18.4, 5.1],
    "selected_atom_count": 10,
    "selected_object_count": 1,
    "chain": "A",
    "objects": ["residue:A:57:HIS"],
    "secondary_structure": "alpha_helix",
    "exposure_status": "partially_buried",
    "interface_chains": ["B"],
    "structural_profile": {
      "secondary_structure": "alpha_helix",
      "secondary_structure_code": "H",
      "secondary_structure_source": "mkdssp",
      "exposure": "partially_buried",
      "exposure_source": "freesasa",
      "surface_status": "intermediate",
      "sasa_A2": 42.7,
      "relative_sasa": 0.18,
      "local_contact_count": 15,
      "local_packing": "medium_local_packing",
      "local_packing_source": "local_packing_density",
      "interface_chains": ["B"],
      "ligand_contacts": ["ligand:A:401:ATP"],
      "ligand_contact_shell_ids": ["ligand_contact_shell:A:401:ATP"],
      "pocket_ids": [],
      "interchain_contact_interface_ids": ["interchain_contact_interface:A-B"],
      "biological_interface_ids": ["biological_interface:pisa:1"],
      "secondary_structure_element": "secondary_structure:A:alpha_helix:52-65"
    },
    "near_ligands": [
      {
        "object_id": "ligand:A:401:ATP",
        "min_distance_A": 6.8,
        "annotations": {
          "local_packing": "low_local_packing",
          "local_packing_source": "local_packing_density"
        }
      }
    ],
    "distance_to_chain_centroid_A": 11.2,
    "plain_language": "HIS57 is in chain A. Secondary structure is alpha_helix. Exposure is partially_buried (freesasa)."
  },
  "evidence": [
    {"type": "metric", "metric": "min_distance_to_ATP", "value": 6.8, "unit": "angstrom"}
  ],
  "limitations": [
    "Interface status uses inter-chain heavy-atom contacts, not biological assembly metadata."
  ]
}
```

`exposure_source="freesasa"` means exact SASA was computed by the FreeSASA
backend. `local_packing_source="local_packing_density"` is a residue-center
packing signal; it is not solvent exposure and includes a limitation when exact
SASA is unavailable.

---

## 5. `context()`

Retrieve compact local/global context relevant to a selection.

```python
session.context(
    selection: str,
    frame: int | str = 0,
    scale: str | None = None,
    radius: float | None = None,
    budget: str | None = None,
    focus: str | list[str] | None = None,
)
```

`scale` is a transparent preset for common inspection depths. It expands to
ordinary `radius`, `budget`, and `focus` values; explicit `radius`, `budget`, or
`focus` arguments override the preset. The effective controls are returned in the
`ContextResult`.

Supported context scales:

| Scale | Expands to | Use |
| --- | --- | --- |
| `chemical_contacts` | `radius=4.0`, `budget="small"`, `focus="contacts"` | Close contacts around a residue, ligand, atom, or site. |
| `residue_environment` | `radius=8.0`, `budget="medium"`, `focus="general"` | Broader local environment around a residue or region. |
| `ligand_binding_site` | `radius=4.0`, `budget="medium"`, `focus="ligand_contact_shell"` | Ligand/ion contact shell and binding-site residues. |
| `metal_coordination` | `radius=3.0`, `budget="small"`, `focus="metal_coordination"` | Metal coordination around a metal atom or cofactor. |
| `protein_interface` | `radius=5.0`, `budget="medium"`, `focus="interchain_interfaces"` | Inter-chain contacts around a protein chain or interface region. |
| `broad_environment` | `radius=12.0`, `budget="large"`, `focus="general"` | Wider orientation pass before narrowing the inspection. |

Examples:

```python
session.context("A:87", scale="chemical_contacts")
session.context("resname HEM", scale="ligand_binding_site")
session.context("name FE", scale="metal_coordination")
session.context("protein and chain A", scale="protein_interface")
session.context("A:87", scale="residue_environment", focus="hydrogen_bonds")
```

For protein-interface inspection, prefer `protein and chain A` over bare
`chain A` when a PDB file assigns the same chain ID to waters, ligands, or ions.

Supported focus values:

| Focus | Meaning |
| --- | --- |
| `general` | Default distance and relation ranking. |
| `contacts` | Prioritize non-`near` relations. |
| `ligand_contact_shell` | Prioritize ligands and ions. |
| `metal_coordination` | Prioritize `metal_coordination` relations. |
| `interchain_interfaces` | Prioritize inter-chain neighbors. |
| `hydrogen_bonds` | Prioritize `hydrogen_bond` and `polar_contact_candidate` relations. |
| `salt_bridges` | Prioritize `salt_bridge` relations. |
| `hydrophobic_contacts` | Prioritize `hydrophobic_contact` relations. |
| `pi_stacking` | Prioritize `pi_stacking` relations. |
| `water_bridges` | Prioritize `water_bridge_candidate` relations. |
| `steric_clashes` | Prioritize `steric_clash` relations. |

`focus` is an explicit ranking control, not natural-language input. Unknown focus
values are rejected.

Context scale presets are also discoverable from Python:

```python
from molinspect import context_scales

context_scales()
session.context_scales()
```

Both return `ContextScalesResult`:

```json
{
  "count": 6,
  "scales": [
    {
      "id": "ligand_binding_site",
      "radius_A": 4.0,
      "budget": "medium",
      "focus": ["ligand_contact_shell"],
      "description": "Ligand/ion contact shell and nearby binding-site residues."
    }
  ],
  "limitations": [
    "Context scales are retrieval presets, not structural-biology definitions."
  ]
}
```

Timeline metric controls are discoverable the same way:

```python
from molinspect import timeline_metrics

timeline_metrics()
session.timeline_metrics()
```

Both return `TimelineMetricsResult` with each metric ID, required arguments,
example call, expected summary keys, sampled-value keys, event types, and
limitations.

Supported budgets:

```text
small   25 nearby objects
medium  75 nearby objects
large   200 nearby objects
```

### Return

```json
{
  "selection": "chain A and resid 57",
  "frame": 0,
  "scale": "residue_environment",
  "radius_A": 8.0,
  "budget": "medium",
  "focus": ["hydrogen_bonds"],
  "objects": [
    {"id": "residue:A:57:HIS", "type": "residue", "name": "HIS57"},
    {
      "id": "residue:A:102:ASP",
      "type": "residue",
      "name": "ASP102",
      "annotations": {"secondary_structure": "loop", "exposure": "exposed"}
    },
    {"id": "ligand:A:401:ATP", "type": "ligand", "name": "ATP"}
  ],
  "relations": [
    {
      "source": "residue:A:57:HIS",
      "target": "residue:A:102:ASP",
      "type": "hydrogen_bond",
      "category": "polar",
      "confidence": "geometry",
      "backend": "molinspect_heuristic",
      "definition_id": "hydrogen_bond",
      "definition_source": "literature_informed_heuristic",
      "reference_keys": ["plip_docs", "arpeggio_paper"],
      "min_distance_A": 3.2,
      "cutoff_A": 4.1,
      "angle_deg": 161.4,
      "source_atom": "A:57:HIS:ND1",
      "target_atom": "A:102:ASP:OD1"
    },
    {
      "source": "residue:A:57:HIS",
      "target": "ligand:A:401:ATP",
      "type": "near",
      "min_distance_A": 6.8
    }
  ],
  "summary": "2 nearby objects within 8 A; 1 is a contact at <= 4 A. Relation types: 1 hydrogen_bond, 1 near."
}
```

---

## 6. `timeline()`

Summarize a metric or relation over time.

```python
session.timeline(
    metric: str,
    selection: str | None = None,
    selection1: str | None = None,
    selection2: str | None = None,
    frames: str | tuple[int | str, int | str] = "all",
    stride: int | None = None,
)
```

Supported metrics:

| Metric | Meaning |
| --- | --- |
| `distance` | Minimum heavy-atom distance over frames. |
| `contact` | Boolean contact state over frames using the 4.0 A cutoff. |
| `relation` | Closest-pair relation type and atom-pair evidence over frames; explicit-hydrogen hydrogen bonds use MDAnalysis HydrogenBondAnalysis evidence when available. |
| `hydrogen_bonds` | Explicit-hydrogen donor-acceptor hydrogen-bond occupancy over frames, with donor, hydrogen, acceptor, distance, angle, event-frame, and limitation evidence. |
| `interaction_persistence` | Combined contact, relation, and explicit-H-bond persistence between two selections. |
| `ligand_stability` | Site-aligned ligand RMSD plus contact/H-bond persistence and a structural stability class. |
| `states` | Simple representative conformational states from aligned RMSD and selection-spread features. |
| `centroid_distance` | Center-of-geometry distance between two selections over frames; useful for domain, loop, partner, or site opening/closing motions. |
| `rmsd` | Kabsch-aligned RMSD for one `selection` against the first selected frame. |
| `rmsf` | Kabsch-aligned per-atom RMSF summarized by residue/object for one `selection`. |
| `mobility` | Mobility-oriented RMSF summary for one `selection`. |
| `displacement` | Kabsch-aligned per-residue/object displacement from the first selected frame. |
| `selection_spread` | Radius-of-gyration spread proxy for selected atoms. |

For programmatic discovery, prefer:

```python
from molinspect import timeline_metrics

guide = timeline_metrics()
guide.metrics[0].required_arguments
guide.metrics[0].summary_keys
```

`frames` can be `"all"` or an inclusive `(start, end)` tuple, for example
`frames=(0, -1)` or `frames=("first", "last")`.

### Example

```python
session.timeline(metric="distance", selection1="chain A and resid 57", selection2="resname ATP")
session.timeline(metric="relation", selection1="chain A and resname HEM and name FE", selection2="chain A and resid 87")
session.timeline(metric="hydrogen_bonds", selection1="chain A and resid 57", selection2="resname ATP")
session.timeline(metric="interaction_persistence", selection1="chain A and resid 57", selection2="resname ATP")
session.timeline(metric="ligand_stability", selection1="resname ATP", selection2="around 4 of resname ATP and protein")
session.timeline(metric="states", selection="chain A")
session.timeline(metric="centroid_distance", selection1="resid 1-29", selection2="resid 122-159")
session.timeline(metric="rmsd", selection="chain A")
session.timeline(metric="mobility", selection="ligand_contact_shell:A:142:HEM")
session.timeline(metric="selection_spread", selection="ligand_contact_shell:A:142:HEM")
```

### Return

```json
{
  "metric": "distance",
  "selection1": "chain A and resid 57",
  "selection2": "resname ATP",
  "frames_analyzed": 1000,
  "summary_text": "distance over 1000 frames: min 3.1 A, median 6.4 A, max 14.2 A; contact occupancy 0.43; dominant relation near.",
  "summary": {
    "min_A": 3.1,
    "median_A": 6.4,
    "max_A": 14.2,
    "contact_occupancy": 0.43,
    "relation_type_counts": {"hydrogen_bond": 400, "near": 600},
    "relation_occupancy": {"hydrogen_bond": 0.4, "near": 0.6},
    "dominant_relation_type": "near",
    "representative_frames": {
      "closest": {"frame": 120, "distance_A": 3.1, "relation_type": "hydrogen_bond"},
      "farthest": {"frame": 450, "distance_A": 14.2, "relation_type": "near"}
    }
  },
  "events": [
    {"type": "contact_forms", "frame": 120, "distance_A": 3.8},
    {"type": "contact_breaks", "frame": 450, "distance_A": 5.1}
  ],
  "sampled_values": [
    {"frame": 0, "value_A": 6.8},
    {"frame": 100, "value_A": 4.1}
  ]
}
```

For `metric="relation"`, `sampled_values` include closest atom-pair evidence:

```json
[
  {
    "frame": 0,
    "relation_type": "metal_coordination",
    "distance_A": 2.143,
    "source_atom": "A:142:HEM:FE",
    "target_atom": "A:87:HIS:NE2"
  }
]
```

For `metric="hydrogen_bonds"`, the summary reports occupancy and top
donor-acceptor pairs:

```json
{
  "backend": "MDAnalysis HydrogenBondAnalysis",
  "hydrogen_bond_occupancy": 0.5,
  "total_hydrogen_bond_observations": 1,
  "donor_acceptor_distance_cutoff_A": 3.0,
  "donor_hydrogen_acceptor_angle_min_deg": 150.0,
  "top_hydrogen_bond_pairs": [
    {
      "donor_atom": "A:57:SER:OG",
      "acceptor_atom": "A:101:ATP:O1",
      "frame_count": 50,
      "occupancy": 0.5
    }
  ]
}
```

This metric requires explicit hydrogens or donor-hydrogen topology that
MDAnalysis can resolve. When that evidence is absent, the result returns zero
observations and a limitation instead of inferring hydrogen bonds from heavy
atoms alone.

For `metric="interaction_persistence"`, the summary combines the evidence users
usually need for a residue-pair or ligand-site relation:

```json
{
  "contact_occupancy": 0.72,
  "hydrogen_bond_occupancy": 0.31,
  "dominant_relation_type": "hydrogen_bond",
  "relation_occupancy": {"hydrogen_bond": 0.31, "near": 0.69},
  "min_distance_A": 2.8,
  "max_distance_A": 8.1,
  "top_hydrogen_bond_pairs": []
}
```

For `metric="ligand_stability"`, `selection1` should be the ligand and
`selection2` should be the binding-site atoms or local site selection. The
summary includes site-aligned ligand RMSD, contact/H-bond occupancy, and
`stability_class`. This class is a structural summary rule, not a binding-free
energy estimate.

For `metric="states"`, `selection` is clustered into simple representative frame
states using aligned RMSD and selection-spread features. The output is intended
for inspection and representative-frame choice, not kinetic state modeling.

---

## 7. `compare()`

Compare a selection across two frames/states.

```python
session.compare(
    selection: str,
    frame_a: int | str,
    frame_b: int | str,
    radius: float = 8.0,
)
```

### Return

```json
{
  "selection": "around 8 of resname ATP",
  "frame_a": 0,
  "frame_b": 999,
  "main_changes": [
    {
      "type": "selection_rmsd",
      "from_frame": 0,
      "to_frame": 999,
      "rmsd_A": 2.4,
      "alignment": "kabsch_over_selection"
    },
    {
      "type": "distance_change",
      "object": "ligand:A:401:ATP",
      "from_A": 3.8,
      "to_A": 9.2,
      "delta_A": 5.4
    },
    {
      "type": "contact_lost",
      "object": "ligand:A:401:ATP",
      "from_relation": "nonbonded_contact",
      "to_relation": "near"
    }
  ],
  "summary": "0 nearby objects gained, 0 lost, 4 shared within 8 A; selection RMSD 2.4 A."
}
```
