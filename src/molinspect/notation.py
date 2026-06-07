"""Public notation help shared by API boundary errors and docs.

These strings describe the stable MolInspect-facing subset. The resolver may
delegate to MDAnalysis internally, but user-facing calls should stay within this
small, documented notation.
"""

SELECTION_NOTATION_HELP = (
    "Use explicit selection notation such as 'chain A and resid 57', 'A:57', "
    "'resid 10-30', 'resname ATP', 'name FE', 'ligand', 'water', 'ion', or "
    "'around 4 of chain A and resid 57'. Combine selectors with 'and', 'or', "
    "'not', and parentheses."
)

FRAME_NOTATION_HELP = (
    "Use a zero-based integer frame index, a negative index from the end, "
    "or the named frames 'first' and 'last'."
)

FRAME_RANGE_NOTATION_HELP = (
    "Use frames='all' or an inclusive (start, end) tuple where each endpoint "
    "uses frame notation."
)

OBJECT_TYPE_NOTATION_HELP = (
    "Supported object types are 'atom', 'residue', 'ligand', 'ion', 'water', "
    "'chain', 'secondary_structure', 'loop', 'ligand_contact_shell', 'pocket', "
    "'interchain_contact_interface', and 'biological_interface'. Plural aliases such "
    "as 'residues' are accepted."
)

CONTAINS_NOTATION_HELP = (
    "contains must be a non-empty literal substring such as 'ATP', 'HEM', "
    "'A:57', or 'ligand:A:142:HEM'. It is not natural-language search."
)

CONTEXT_FOCUS_NOTATION_HELP = (
    "Supported focus values are 'general', 'contacts', 'ligand_contact_shell', "
    "'metal_coordination', 'interchain_interfaces', 'hydrogen_bonds', 'salt_bridges', "
    "'hydrophobic_contacts', 'pi_stacking', 'water_bridges', and 'steric_clashes'."
)

CONTEXT_SCALE_NOTATION_HELP = (
    "Supported context scale values are 'chemical_contacts', 'residue_environment', "
    "'ligand_binding_site', 'metal_coordination', 'protein_interface', and "
    "'broad_environment'. A scale expands to explicit radius, budget, and focus values."
)

BUDGET_NOTATION_HELP = "budget must be one of 'small', 'medium', or 'large'."

TIMELINE_METRIC_NOTATION_HELP = (
    "timeline metric must be one of 'distance', 'contact', 'relation', "
    "'hydrogen_bonds', 'interaction_persistence', 'ligand_stability', 'states', "
    "'centroid_distance', 'rmsd', 'rmsf', 'mobility', 'displacement', "
    "or 'selection_spread'."
)
