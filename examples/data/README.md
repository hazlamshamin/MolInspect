# Real PDB Example Data

These structures are small, well-known examples for exercising MolInspect's first API
surface with real molecular files.

| PDB ID | Why it is useful | Source |
|---|---|---|
| `1CRN` | Tiny single-chain crambin structure; good for simple loading and residue selections. | https://files.rcsb.org/download/1CRN.pdb |
| `1UBQ` | Canonical compact ubiquitin structure; good for a familiar single-chain protein. | https://files.rcsb.org/download/1UBQ.pdb |
| `4HHB` | Classic hemoglobin tetramer with heme ligands; good for chain and ligand object checks. | https://files.rcsb.org/download/4HHB.pdb |

Fetch or refresh these files with:

```bash
uv run python scripts/fetch_example_pdbs.py
```

The automated unit tests use generated tiny fixtures by default. The real PDB files are
kept here for local validation, examples, and visual QA.
