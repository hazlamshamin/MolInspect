# API Notation Guide

MolInspect API inputs use explicit structural notation, not natural-language
commands. This page defines the string values that appear in public calls and
outputs.

## Object IDs

Object IDs are stable labels for molecular objects in the loaded target.

```text
atom:{chain}:{resid}{icode}:{resname}:{atom_name}:{atom_index0}
residue:{chain}:{resid}{icode}:{resname}
ligand:{chain}:{resid}{icode}:{resname}
ion:{chain}:{resid}{icode}:{resname}
water:{chain}:{resid}{icode}:{resname}
chain:{chain}
selection_region:{n_residues}res:{first_atom_index0}-{last_atom_index0}:{fingerprint}
secondary_structure:{chain}:{kind}:{start_resid}-{end_resid}
loop:{chain}:{start_resid}-{end_resid}
ligand_contact_shell:{chain}:{resid}{icode}:{resname}
pocket:p2rank:{rank}
pocket:fpocket:{rank}
interchain_contact_interface:{chain_a}-{chain_b}
biological_interface:pisa:{interface_id}
```

Examples:

```text
residue:A:87:HIS
ligand:A:142:HEM
ion:A:301:NA
atom:A:142:HEM:FE:1064
chain:A
ligand_contact_shell:A:142:HEM
pocket:p2rank:1
interchain_contact_interface:A-B
biological_interface:pisa:1
secondary_structure:A:alpha_helix:21-35
```

Field meanings:

- `chain` is the PDB chain ID when present; otherwise MolInspect uses the segment ID, then `_`.
- `resid` is the residue number from the topology. Insertion codes are appended, such as `57A`.
- `resname` is the residue or ligand name, such as `HIS`, `ATP`, or `HEM`.
- `atom_index0` is the zero-based atom index from the loaded topology.
- `icode` is an insertion code appended to `resid` when present, such as `57A`.
- Atom objects also expose `altloc` when the topology provides an alternate-location
  identifier. The atom index remains part of the object ID for uniqueness.
- `selection_region` is a compact synthetic object for broad selections with more than 25 residues.
- `ligand_contact_shell` is a direct ligand/ion heavy-atom contact shell.
- `pocket` is reserved for a successful backend pocket call. P2Rank is preferred
  for ligand-binding-site prediction; fpocket is a geometric fallback. It is not
  a synonym for a ligand contact shell.
- `interchain_contact_interface` is named by the two chains with inter-chain
  heavy-atom contacts.
- `biological_interface` is reserved for successful PISA interface results.
  It carries PISA buried-surface/energy evidence and is not a synonym for
  `interchain_contact_interface`.
- `secondary_structure` and `loop` IDs are emitted when DSSP annotations are available.

## Atom Labels

Relation evidence uses compact atom labels:

```text
{chain}:{resid}{icode}:{resname}:{atom_name}
```

Example:

```text
A:87:HIS:NE2
```

Atom labels are evidence labels, not stable object IDs. Use object IDs when
referencing atoms as objects.

## Selection Notation

Selections identify what the inspection should operate on. The stable public
subset is:

```text
chain A
resid 87
resid 10-30
resname HEM
name FE
protein
nucleic
ligand
water
ion
A:87
around 4 of chain A and resid 87
```

Object IDs returned by `objects()`, `locate()`, or `context()` can also be passed
back as selections when they identify a concrete object or backend landmark:

```text
residue:A:87:HIS
ligand:A:142:HEM
atom:A:142:HEM:FE:1064
ligand_contact_shell:A:142:HEM
pocket:p2rank:1
biological_interface:pisa:1
```

Selectors can be combined with:

```text
and
or
not
(...)
```

Examples:

```python
session.locate("chain A and resid 87")
session.context("chain A and resname HEM and name FE", radius=3.0)
session.timeline(metric="relation", selection1="A:87", selection2="resname HEM")
session.timeline(metric="rmsd", selection="chain A")
session.timeline(metric="mobility", selection="chain A")
session.resolve_selection("around 4 of chain A and resid 87")
session.resolve_selection("ligand_contact_shell:A:142:HEM")
session.resolve_selection("biological_interface:pisa:1")
```

`A:87` is shorthand for `chain A and resid 87`. `around 4 of ...` selects atoms
within 4 Angstrom of the inner selection.

## Frames

Frame inputs use zero-based indexing:

```text
0       first frame
1       second frame
-1      last frame
"first" first frame
"last"  last frame
```

Timeline frame ranges use:

```python
frames="all"
frames=(0, -1)
frames=("first", "last")
```

Frame range tuples are inclusive. Unsupported frame labels raise `FrameError`.

## Units And Field Suffixes

- `_A` means Angstrom, for example `min_distance_A`.
- `occupancy` and `contact_occupancy` are fractions from `0.0` to `1.0`.
- `frame`, `frame_a`, and `frame_b` in outputs are normalized zero-based frame indices.

## Context Scales

`scale` is a transparent preset for common inspection depth. It does not parse
language. It only chooses concrete `radius`, `budget`, and `focus` values.

```text
chemical_contacts    radius 4.0 A, budget small,  focus contacts
residue_environment  radius 8.0 A, budget medium, focus general
ligand_binding_site  radius 4.0 A, budget medium, focus ligand_contact_shell
metal_coordination   radius 3.0 A, budget small,  focus metal_coordination
protein_interface    radius 5.0 A, budget medium, focus interchain_interfaces
broad_environment    radius 12.0 A, budget large, focus general
```

Examples:

```python
session.context("A:87", scale="chemical_contacts")
session.context("resname HEM", scale="ligand_binding_site")
session.context("name FE", scale="metal_coordination")
session.context("protein and chain A", scale="protein_interface")
```

Explicit controls override the preset while the returned result still reports
the selected scale and the effective controls:

```python
session.context("A:87", scale="residue_environment", radius=5.0, focus="hydrogen_bonds")
```

For protein-interface inspection, prefer `protein and chain A` over bare
`chain A` when a PDB file assigns the same chain ID to waters, ligands, or ions.

The same list is available programmatically:

```python
from molinspect import context_scales

context_scales()
session.context_scales()
```

Timeline metric call shapes are also available programmatically:

```python
from molinspect import timeline_metrics

timeline_metrics()
session.timeline_metrics()
```

## Context Focus

`focus` ranks the returned context evidence. It does not parse language.

```text
general              default ranking
contacts             prioritize non-near relations
ligand_contact_shell prioritize ligands and ions
metal_coordination   prioritize metal_coordination relations
interchain_interfaces prioritize inter-chain neighbors
hydrogen_bonds       prioritize hydrogen_bond and polar_contact_candidate relations
salt_bridges         prioritize salt_bridge relations
hydrophobic_contacts prioritize hydrophobic_contact relations
pi_stacking          prioritize pi_stacking relations
water_bridges        prioritize water_bridge_candidate relations
steric_clashes       prioritize steric_clash relations
```

## Budgets

`budget` caps how many nearby objects `context()` returns:

```text
small    25 objects
medium   75 objects
large    200 objects
```

## Timeline Metrics

```text
distance   minimum heavy-atom distance over frames
contact    boolean contact state over frames, using the 4.0 A cutoff
relation   closest-pair relation type and atom-pair evidence over frames;
           explicit-hydrogen hydrogen bonds use MDAnalysis backend evidence
hydrogen_bonds explicit-hydrogen donor-acceptor H-bond occupancy over frames
interaction_persistence combined contact, relation, and H-bond persistence summary
ligand_stability site-aligned ligand RMSD plus contact/H-bond persistence summary
states     simple representative conformational states from RMSD/spread features
centroid_distance distance between centers of geometry of two selections over frames
rmsd       Kabsch-aligned RMSD of one selection against the first selected frame
rmsf       Kabsch-aligned per-atom RMSF summarized by residue/object
mobility   alias-style mobility view over RMSF summaries
displacement Kabsch-aligned per-residue/object displacement from the first selected frame
selection_spread radius-of-gyration spread proxy for selected atoms
```

## Relation Types

Relation labels describe current structural evidence. Outputs carry attribution
fields such as `definition_id`, `definition_source`, `backend`,
`reference_keys`, `cutoff_A`, and `limitations`. The full registry lives in
`src/molinspect/definitions/registry.py` for users who need implementation-level
details.

```text
topology_bond           covalent bond read from topology/connectivity records
inferred_covalent_bond  local covalent-bond estimate when topology lacks bonds
steric_clash            non-bonded heavy atoms with vdW overlap above the clash threshold
metal_coordination      metal/coordinating atom pair within 3.0 A
salt_bridge             oppositely charged atoms within 5.5 A
hydrogen_bond           donor/acceptor pair within 4.1 A with explicit H angle >= 100 degrees
polar_contact_candidate polar atom pair within 4.1 A without donor/acceptor-angle validation
pi_stacking             aromatic ring pair passing center, plane-angle, and offset checks
hydrophobic_contact     hydrophobic side-chain carbon pair within 4.0 A
water_bridge_candidate  polar atoms bridged by one water oxygen with 2.5-4.1 A legs
nonbonded_contact       closest non-covalent heavy-atom distance <= 4.0 A
near                    within the requested context radius, but not a contact
```

Relation outputs include `category`, `confidence`, `backend`, `definition_id`,
`definition_source`, `reference_keys`, `cutoff_A`, optional `angle_deg`, atom
labels, metric evidence, and limitations. `_candidate` is reserved for cases
where the geometry lacks a required validation signal, such as missing explicit
hydrogen atoms for hydrogen-bond angle checks.
