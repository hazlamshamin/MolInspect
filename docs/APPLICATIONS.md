# Applications

MolInspect is meant for workflows that need compact, evidence-backed structural
context from existing molecular data.

## Structural Biology

Use MolInspect to turn a local region into a structured evidence card:

- identify nearby residues, ligands, ions, waters, pockets, and interfaces;
- report secondary-structure and exposure context;
- expose atom-pair distances, relation labels, cutoffs, backends, and
  limitations;
- keep output small enough for notes, reports, and LLM context.

Example:

```python
session.context("chain A and resname HEM", scale="ligand_binding_site")
```

## Protein Engineering

Use MolInspect before mutation or design work to understand what a residue is
doing structurally:

- whether it is buried, exposed, ligand-facing, pocket-facing, or interface-facing;
- whether it contributes salt-bridge, metal-coordination, polar, hydrophobic, or
  aromatic evidence;
- which local residues and ligands would be affected by a proposed mutation.

Example:

```python
session.locate("chain A and resid 87")
session.context("chain A and resid 87", scale="chemical_contacts")
```

## Ligand And Cofactor Analysis

Use MolInspect to summarize binding-site evidence without reading the full PDB:

- ligand-contact-shell objects;
- nearby residue keys and object IDs;
- metal coordination and polar-contact candidates;
- pocket objects when a pocket backend is available.

Example:

```python
session.objects(type="ligand_contact_shell", contains="HEM")
session.context("chain A and resname HEM and name FE", radius=3.0,
                focus="metal_coordination")
```

## Protein-Protein Interfaces

Use MolInspect to inspect chain interfaces and biological-interface annotations:

- inter-chain contact interface objects;
- PISA-backed biological interfaces when available;
- interface-aware context slices instead of whole-chain dumps.

Example:

```python
session.objects(type="interchain_contact_interface")
session.context("protein and chain A", scale="protein_interface")
```

## Molecular Dynamics And Ensembles

Use MolInspect to summarize structural change over frames:

- distance/contact timelines;
- relation occupancy;
- explicit-hydrogen hydrogen-bond occupancy with donor, hydrogen, acceptor,
  distance, angle, and event-frame evidence;
- interaction-persistence summaries that combine contact, relation, and H-bond
  evidence;
- ligand/site stability summaries from site-aligned ligand RMSD and interaction
  persistence;
- simple representative conformational-state summaries;
- Kabsch-aligned RMSD and RMSF-style mobility;
- displacement and spread summaries;
- representative frames and event-like changes.

Example:

```python
session.timeline(metric="relation", selection1="chain A and resid 57",
                 selection2="resname ATP")
session.timeline(metric="hydrogen_bonds", selection1="chain A and resid 57",
                 selection2="resname ATP")
session.timeline(metric="interaction_persistence", selection1="chain A and resid 57",
                 selection2="resname ATP")
session.timeline(metric="ligand_stability", selection1="resname ATP",
                 selection2="around 4 of resname ATP and protein")
session.timeline(metric="states", selection="chain A")
session.timeline(metric="mobility", selection="chain A")
session.compare("around 8 of resname ATP", frame_a="first", frame_b="last")
```

## LLM And Agent Workflows

Use MolInspect as a high-signal retrieval layer before reasoning:

- retrieve only the relevant structural slice;
- keep object IDs stable for follow-up calls;
- preserve exact measurements and limitations in typed outputs;
- optionally render examples for human validation.

MolInspect should be treated as the source of structural evidence. The LLM or
agent should reason from returned objects, relations, metrics, frames, and
limitations instead of from raw coordinate dumps.
