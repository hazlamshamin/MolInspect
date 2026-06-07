"""Discoverable timeline metric specifications."""

from __future__ import annotations

from dataclasses import dataclass

from ..schemas import TimelineMetricInfo, TimelineMetricsResult


@dataclass(frozen=True, slots=True)
class TimelineMetricSpec:
    """Concrete call/return guide for one timeline metric."""

    metric: str
    required_arguments: tuple[str, ...]
    best_for: str
    description: str
    summary_keys: tuple[str, ...]
    sampled_value_keys: tuple[str, ...]
    event_types: tuple[str, ...]
    example: str
    limitations: tuple[str, ...] = ()


TIMELINE_METRIC_SPECS: tuple[TimelineMetricSpec, ...] = (
    TimelineMetricSpec(
        metric="distance",
        required_arguments=("selection1", "selection2"),
        best_for="Minimum heavy-atom distance between two selections over frames.",
        description="Use for closest-approach questions and contact-distance evidence.",
        summary_keys=("min_A", "median_A", "max_A", "mean_A", "contact_occupancy"),
        sampled_value_keys=("frame", "value_A"),
        event_types=("contact_forms", "contact_breaks", "ligand_moved_toward", "ligand_moved_away"),
        example=(
            'session.timeline(metric="distance", selection1="chain A and resid 57", '
            'selection2="resname ATP")'
        ),
        limitations=("Distance is a minimum heavy-atom distance, not a center distance.",),
    ),
    TimelineMetricSpec(
        metric="contact",
        required_arguments=("selection1", "selection2"),
        best_for="Boolean contact persistence between two selections.",
        description="Use when contact occupancy and form/break frames matter more than relation type.",
        summary_keys=("contact_occupancy", "contact_cutoff_A", "representative_frames"),
        sampled_value_keys=("frame", "contact"),
        event_types=("contact_forms", "contact_breaks"),
        example=(
            'session.timeline(metric="contact", selection1="chain A and resid 57", '
            'selection2="resname ATP")'
        ),
        limitations=("Contact uses the registered heavy-atom cutoff.",),
    ),
    TimelineMetricSpec(
        metric="relation",
        required_arguments=("selection1", "selection2"),
        best_for="Typed closest relation and atom-pair evidence over frames.",
        description=(
            "Use for relation changes such as metal coordination, salt bridge, "
            "hydrogen bond, steric clash, or near/contact state."
        ),
        summary_keys=("relation_occupancy", "dominant_relation_type", "representative_frames"),
        sampled_value_keys=("frame", "relation_type", "distance_A", "source_atom", "target_atom"),
        event_types=(
            "relation_formed",
            "relation_lost",
            "relation_changed",
            "contact_formed",
            "contact_lost",
            "contact_forms",
            "contact_breaks",
        ),
        example=(
            'session.timeline(metric="relation", selection1="chain A and resname HEM and name FE", '
            'selection2="chain A and resid 87")'
        ),
        limitations=("Non-H-bond relation types use closest-pair structural evidence.",),
    ),
    TimelineMetricSpec(
        metric="hydrogen_bonds",
        required_arguments=("selection1", "selection2"),
        best_for="Explicit-hydrogen donor-acceptor H-bond occupancy.",
        description="Use when the topology has hydrogens and donor/hydrogen/acceptor evidence is needed.",
        summary_keys=(
            "hydrogen_bond_occupancy",
            "total_hydrogen_bond_observations",
            "top_hydrogen_bond_pairs",
        ),
        sampled_value_keys=("frame", "hydrogen_bond_count", "has_hydrogen_bond", "observations"),
        event_types=("hydrogen_bond_forms", "hydrogen_bond_breaks"),
        example=(
            'session.timeline(metric="hydrogen_bonds", selection1="chain A and resid 57", '
            'selection2="resname ATP")'
        ),
        limitations=("Requires explicit hydrogens or donor-hydrogen topology MDAnalysis can resolve.",),
    ),
    TimelineMetricSpec(
        metric="interaction_persistence",
        required_arguments=("selection1", "selection2"),
        best_for="One high-level persistence card for two structural selections.",
        description=(
            "Combines contact occupancy, relation occupancy, explicit-H-bond occupancy, "
            "events, and representative frames."
        ),
        summary_keys=(
            "contact_occupancy",
            "hydrogen_bond_occupancy",
            "dominant_relation_type",
            "relation_occupancy",
            "min_distance_A",
            "max_distance_A",
        ),
        sampled_value_keys=(
            "frame",
            "relation_type",
            "distance_A",
            "contact",
            "hydrogen_bond_count",
        ),
        event_types=("contact_forms", "contact_breaks", "hydrogen_bond_forms", "hydrogen_bond_breaks"),
        example=(
            'session.timeline(metric="interaction_persistence", selection1="chain A and resid 57", '
            'selection2="resname ATP")'
        ),
        limitations=("This is an evidence summary, not a binding-energy estimate.",),
    ),
    TimelineMetricSpec(
        metric="ligand_stability",
        required_arguments=("selection1", "selection2"),
        best_for="Ligand pose stability inside an explicit site selection.",
        description=(
            "Treat selection1 as the ligand and selection2 as the binding site; reports "
            "site-aligned ligand RMSD plus interaction persistence."
        ),
        summary_keys=(
            "stability_class",
            "max_site_aligned_ligand_rmsd_A",
            "contact_occupancy",
            "hydrogen_bond_occupancy",
            "representative_frames",
        ),
        sampled_value_keys=(
            "frame",
            "site_aligned_ligand_rmsd_A",
            "ligand_site_centroid_distance_A",
            "relation_type",
            "hydrogen_bond_count",
        ),
        event_types=("ligand_pose_repositioned", "contact_forms", "contact_breaks"),
        example=(
            'session.timeline(metric="ligand_stability", selection1="resname ATP", '
            'selection2="around 4 of resname ATP and protein")'
        ),
        limitations=("Selection1 must keep consistent ligand atom correspondence across frames.",),
    ),
    TimelineMetricSpec(
        metric="states",
        required_arguments=("selection",),
        best_for="Simple representative conformational states for one selection.",
        description="Bins frames from aligned RMSD and selection-spread features for inspection.",
        summary_keys=("state_count", "state_method", "rmsd_range_A", "states", "representative_frames"),
        sampled_value_keys=("frame", "state_id", "aligned_rmsd_A", "selection_spread_A"),
        event_types=("state_transition",),
        example='session.timeline(metric="states", selection="chain A")',
        limitations=("Deterministic inspection bins, not kinetic clustering or an MSM.",),
    ),
    TimelineMetricSpec(
        metric="centroid_distance",
        required_arguments=("selection1", "selection2"),
        best_for="Opening/closing movement between two regions, domains, loops, or partners.",
        description="Measures center-of-geometry distance over frames.",
        summary_keys=("min_centroid_distance_A", "max_centroid_distance_A", "representative_frames"),
        sampled_value_keys=("frame", "centroid_distance_A"),
        event_types=("centroid_distance_increased", "centroid_distance_decreased"),
        example=(
            'session.timeline(metric="centroid_distance", selection1="resid 1-29", '
            'selection2="resid 122-159")'
        ),
        limitations=("Centroid distance is not a contact distance and is not mass-weighted.",),
    ),
    TimelineMetricSpec(
        metric="rmsd",
        required_arguments=("selection",),
        best_for="Aligned conformational deviation from the first selected frame.",
        description="Kabsch-aligns the selected atoms against the first selected frame.",
        summary_keys=("min_rmsd_A", "median_rmsd_A", "max_rmsd_A", "reference_frame"),
        sampled_value_keys=("frame", "rmsd_A"),
        event_types=("rmsd_peak",),
        example='session.timeline(metric="rmsd", selection="chain A")',
        limitations=("Requires consistent atom correspondence across frames.",),
    ),
    TimelineMetricSpec(
        metric="rmsf",
        required_arguments=("selection",),
        best_for="Per-residue/object mobility from aligned atom RMSF.",
        description="Summarizes Kabsch-aligned per-atom RMSF by selected residue/object.",
        summary_keys=("min_rmsf_A", "median_rmsf_A", "max_rmsf_A", "most_mobile_object"),
        sampled_value_keys=("object", "mean_rmsf_A", "max_atom_rmsf_A", "selected_atom_count"),
        event_types=("mobility_peak",),
        example='session.timeline(metric="rmsf", selection="chain A")',
        limitations=("RMSF requires consistent atom correspondence across frames.",),
    ),
    TimelineMetricSpec(
        metric="mobility",
        required_arguments=("selection",),
        best_for="User-facing mobility summary over residues or objects.",
        description="Same RMSF evidence as rmsf, named for inspection workflows.",
        summary_keys=("min_rmsf_A", "median_rmsf_A", "max_rmsf_A", "most_mobile_object"),
        sampled_value_keys=("object", "mean_rmsf_A", "max_atom_rmsf_A", "selected_atom_count"),
        event_types=("mobility_peak",),
        example='session.timeline(metric="mobility", selection="chain A")',
        limitations=("Mobility is RMSF-like evidence, not a thermodynamic flexibility model.",),
    ),
    TimelineMetricSpec(
        metric="displacement",
        required_arguments=("selection",),
        best_for="Which residues/objects moved most from the first selected frame.",
        description="Reports aligned per-frame displacement summaries and peak frames.",
        summary_keys=("min_displacement_A", "median_displacement_A", "max_displacement_A"),
        sampled_value_keys=(
            "frame",
            "mean_displacement_A",
            "max_displacement_A",
            "most_displaced_object",
        ),
        event_types=("displacement_peak",),
        example='session.timeline(metric="displacement", selection="chain A")',
        limitations=("Per-residue displacement is aggregated from selected atoms only.",),
    ),
    TimelineMetricSpec(
        metric="selection_spread",
        required_arguments=("selection",),
        best_for="Compactness/opening of one selected atom set.",
        description="Radius-of-gyration-style spread proxy over frames.",
        summary_keys=("min_selection_spread_A", "median_selection_spread_A", "max_selection_spread_A"),
        sampled_value_keys=("frame", "selection_spread_A"),
        event_types=("selection_spread_increased", "selection_spread_decreased"),
        example='session.timeline(metric="selection_spread", selection="ligand_contact_shell:A:142:HEM")',
        limitations=("Spread is a geometry proxy, not a pocket-volume calculation.",),
    ),
)


SUPPORTED_TIMELINE_METRICS = frozenset(spec.metric for spec in TIMELINE_METRIC_SPECS)


def timeline_metric_specs() -> tuple[TimelineMetricSpec, ...]:
    """Return the supported timeline metric specifications."""

    return TIMELINE_METRIC_SPECS


def timeline_metrics() -> TimelineMetricsResult:
    """Return user-facing timeline metric call and output guidance."""

    metrics = [
        TimelineMetricInfo(
            metric=spec.metric,
            required_arguments=list(spec.required_arguments),
            optional_arguments=["frames", "stride"],
            best_for=spec.best_for,
            description=spec.description,
            summary_keys=list(spec.summary_keys),
            sampled_value_keys=list(spec.sampled_value_keys),
            event_types=list(spec.event_types),
            example=spec.example,
            limitations=list(spec.limitations),
        )
        for spec in TIMELINE_METRIC_SPECS
    ]
    return TimelineMetricsResult(
        metrics=metrics,
        count=len(metrics),
        limitations=[
            "Timeline metric IDs are explicit controls, not natural-language commands.",
            "Every metric uses MolInspect selection notation and zero-based frame notation.",
        ],
    )
