"""Static spatial validation helpers for external structural benchmarks.

This module deliberately does not define a MolInspect-native benchmark. It gives
the repo a small, typed surface for converting external benchmark labels into
the residue/object notation used by the public APIs, then scoring API evidence.
"""

from __future__ import annotations

from collections.abc import Iterable, Sequence
from math import dist
from pathlib import Path
from typing import Any, Literal

from pydantic import Field
from scipy.spatial import cKDTree

from molinspect.objects import residue_number_with_icode
from molinspect.schemas import ContextResult, MolInspectModel

ATOM3D_PIP_INTERFACE_CUTOFF_A = 6.0

TaskType = Literal[
    "ligand_binding_residue_retrieval",
    "protein_interface_residue_retrieval",
    "pocket_center_hit",
]


class BenchmarkSource(MolInspectModel):
    """External benchmark/source-of-truth candidate for static spatial tasks."""

    id: str
    name: str
    url: str
    task_types: list[TaskType]
    gold_type: str
    fit_for_molinspect: str
    caveats: list[str] = Field(default_factory=list)
    priority: Literal["primary", "secondary", "stretch"]


class ResidueKey(MolInspectModel):
    """Residue identity key that does not require a three-letter residue name."""

    chain: str
    resid: str
    resname: str | None = None

    @property
    def label(self) -> str:
        """Return the chain/residue key used for scoring."""

        return f"{self.chain}:{self.resid}"


class ResidueRetrievalScore(MolInspectModel):
    """Precision/recall score for residue-set retrieval tasks."""

    true_positives: list[str]
    false_positives: list[str]
    false_negatives: list[str]
    precision: float
    recall: float
    f1: float


class PocketCenterScore(MolInspectModel):
    """Distance-threshold score for a predicted pocket center."""

    distance_A: float
    cutoff_A: float
    hit: bool


class StaticSpatialTask(MolInspectModel):
    """One externally sourced static spatial benchmark task after conversion."""

    task_id: str
    benchmark_id: str
    task_type: TaskType
    structure: str | Path
    selection: str
    gold_residues: list[str] = Field(default_factory=list)
    gold_center_A: tuple[float, float, float] | None = None
    radius_A: float = 4.0
    focus: str | list[str] | None = None
    notes: str | None = None


class BioLiPRecord(MolInspectModel):
    """Parsed fields needed from a BioLiP/BioLiP2 ligand-binding-site row."""

    pdb_id: str
    receptor_chain: str
    resolution_A: float | None = None
    binding_site_code: str | None = None
    ligand_id: str
    ligand_chain: str
    ligand_serial: str
    ligand_auth_resid: str | None = None
    binding_site_residues: list[ResidueKey]

    @property
    def ligand_selection(self) -> str:
        """Return a MolInspect selection for the ligand record."""

        resid = self.ligand_auth_resid or self.ligand_serial
        parts = [f"resid {resid}", f"resname {self.ligand_id}"]
        if self.ligand_chain:
            parts.insert(0, f"chain {self.ligand_chain}")
        return " and ".join(parts)


def external_static_benchmark_sources() -> list[BenchmarkSource]:
    """Return the external source portfolio for static spatial validation."""

    return [
        BenchmarkSource(
            id="biolip2",
            name="BioLiP / BioLiP2 ligand-protein interaction database",
            url="https://zhanggroup.org/BioLiP/",
            task_types=["ligand_binding_residue_retrieval"],
            gold_type="curated biologically relevant ligand-binding residues",
            fit_for_molinspect=(
                "Primary source for checking whether context(ligand, focus='ligand_contact_shell') "
                "retrieves the externally annotated binding-site residues."
            ),
            caveats=[
                "Rows must be filtered to biologically relevant ligands; crystallization additives "
                "and ambiguous ligands should not be used as validation-gate tasks.",
                "BioLiP residue tokens usually identify residue number and amino-acid code, not the "
                "full MolInspect object ID, so scoring uses chain:resid keys.",
            ],
            priority="primary",
        ),
        BenchmarkSource(
            id="cameo_ligand_binding_site",
            name="CAMEO ligand-binding-site assessment protocol",
            url="https://archive.cameo3d.org/cameong_help/3d/",
            task_types=["ligand_binding_residue_retrieval"],
            gold_type="binding-site residue protocol and distance metrics",
            fit_for_molinspect=(
                "Metric reference for deciding whether predicted binding-site residues and "
                "coordinates are close enough to target residues/ligands."
            ),
            caveats=[
                "CAMEO is best used as an assessment protocol/source of target cases, not as the "
                "only gold set.",
                "The protocol is static; it must not be reused as evidence for temporal behavior.",
            ],
            priority="primary",
        ),
        BenchmarkSource(
            id="p2rank_datasets",
            name="P2Rank binding-site benchmark datasets",
            url="https://github.com/rdk/p2rank-datasets",
            task_types=["pocket_center_hit", "ligand_binding_residue_retrieval"],
            gold_type="holo protein-ligand structures used in pocket-prediction benchmarks",
            fit_for_molinspect=(
                "Primary pocket-ranking benchmark family for objects(type='pocket') and pocket "
                "center/lining evidence."
            ),
            caveats=[
                "When MolInspect uses P2Rank as a backend, these datasets validate integration and "
                "evidence presentation more than independent P2Rank accuracy.",
                "Use BioLiP or CAMEO-derived cases alongside this source to avoid backend self-scoring.",
            ],
            priority="primary",
        ),
        BenchmarkSource(
            id="atom3d_pip",
            name="ATOM3D protein interface prediction task",
            url="https://github.com/drorlab/atom3d",
            task_types=["protein_interface_residue_retrieval"],
            gold_type="protein-interface labels from structural complexes",
            fit_for_molinspect=(
                "Primary ML-benchmark source for checking interchain interface residue retrieval "
                "from context(..., focus='interchain_interfaces')."
            ),
            caveats=[
                "ATOM3D packaging is ML-oriented; conversion should keep split identity and avoid "
                "mixing train/test cases.",
                "Interface labels are residue-level, so relation atom-pair evidence should be scored "
                "as secondary evidence.",
            ],
            priority="primary",
        ),
        BenchmarkSource(
            id="docking_benchmark_5",
            name="Protein-Protein Docking Benchmark 5.x",
            url="https://zlab.umassmed.edu/benchmark/",
            task_types=["protein_interface_residue_retrieval"],
            gold_type="curated protein-protein complex benchmark structures",
            fit_for_molinspect=(
                "Independent complex set for biological-interface and interchain-contact retrieval."
            ),
            caveats=[
                "The benchmark is docking-focused; MolInspect should score static interface evidence, "
                "not docking success.",
                "Some cases need biological-assembly handling before judging PISA/interface results.",
            ],
            priority="secondary",
        ),
        BenchmarkSource(
            id="plinder",
            name="PLINDER protein-ligand interaction dataset",
            url="https://plinder.sh/",
            task_types=["ligand_binding_residue_retrieval", "pocket_center_hit"],
            gold_type="large protein-ligand interaction benchmark with modern splits",
            fit_for_molinspect=(
                "Stretch source for larger-scale protein-ligand stress testing after the small "
                "BioLiP/CAMEO/P2Rank gate is reliable."
            ),
            caveats=[
                "Large dataset; avoid making it a default local test.",
                "More docking/ML-oriented than explanation-oriented, so validation gates need narrow "
                "static evidence metrics.",
            ],
            priority="stretch",
        ),
    ]


def residue_key(value: str | ResidueKey) -> ResidueKey:
    """Normalize a residue object ID or chain:resid label into a scoring key."""

    if isinstance(value, ResidueKey):
        return value
    raw = value.strip()
    if raw.startswith("residue:"):
        parts = raw.split(":")
        if len(parts) < 3:
            raise ValueError(f"invalid residue object ID: {value!r}")
        return ResidueKey(
            chain=parts[1],
            resid=parts[2],
            resname=parts[3] if len(parts) > 3 and parts[3] else None,
        )
    parts = raw.split(":")
    if len(parts) in {2, 3}:
        return ResidueKey(
            chain=parts[0],
            resid=parts[1],
            resname=parts[2] if len(parts) == 3 and parts[2] else None,
        )
    raise ValueError(f"expected residue object ID or chain:resid label, got {value!r}")


def context_residue_keys(context: ContextResult) -> list[ResidueKey]:
    """Extract residue keys from a context result's returned objects and relations."""

    labels: dict[str, ResidueKey] = {}
    for obj in context.objects:
        if obj.type == "residue":
            key = residue_key(obj.id)
            labels[key.label] = key
    for relation in context.relations:
        for endpoint in (relation.source, relation.target):
            if endpoint.startswith("residue:"):
                key = residue_key(endpoint)
                labels[key.label] = key
    return list(labels.values())


def score_residue_retrieval(
    predicted: Iterable[str | ResidueKey],
    gold: Iterable[str | ResidueKey],
) -> ResidueRetrievalScore:
    """Score predicted residue identities against an external residue gold set."""

    predicted_labels = {residue_key(value).label for value in predicted}
    gold_labels = {residue_key(value).label for value in gold}
    true_positives = predicted_labels & gold_labels
    false_positives = predicted_labels - gold_labels
    false_negatives = gold_labels - predicted_labels
    precision = _safe_ratio(len(true_positives), len(predicted_labels))
    recall = _safe_ratio(len(true_positives), len(gold_labels))
    f1 = _safe_ratio(2 * precision * recall, precision + recall)
    return ResidueRetrievalScore(
        true_positives=sorted(true_positives),
        false_positives=sorted(false_positives),
        false_negatives=sorted(false_negatives),
        precision=precision,
        recall=recall,
        f1=f1,
    )


def score_pocket_center_hit(
    predicted_center_A: Sequence[float],
    gold_center_A: Sequence[float],
    cutoff_A: float = 4.0,
) -> PocketCenterScore:
    """Score pocket-center distance against a benchmark hit cutoff."""

    if len(predicted_center_A) != 3 or len(gold_center_A) != 3:
        raise ValueError("pocket centers must contain exactly three coordinates")
    if cutoff_A <= 0:
        raise ValueError("cutoff_A must be > 0")
    distance_A = dist(tuple(float(x) for x in predicted_center_A), tuple(float(x) for x in gold_center_A))
    return PocketCenterScore(
        distance_A=round(distance_A, 3),
        cutoff_A=float(cutoff_A),
        hit=distance_A <= cutoff_A,
    )


def atom3d_pip_interface_residues(
    universe: Any,
    source_chain: str,
    target_chain: str,
    cutoff_A: float = ATOM3D_PIP_INTERFACE_CUTOFF_A,
) -> list[ResidueKey]:
    """Return target-chain protein residues using the ATOM3D PPI/PIP interface rule.

    ATOM3D's PPI/PIP README defines interacting amino acids as residues spanning
    two proteins where any heavy atoms are within 6 Angstroms. This helper
    returns the target-chain side of that directed interface.
    """

    if cutoff_A <= 0:
        raise ValueError("cutoff_A must be > 0")
    source = _protein_chain_atomgroup(universe, source_chain)
    target = _protein_chain_atomgroup(universe, target_chain)
    source_heavy = source.select_atoms("not name H*")
    if len(source_heavy) == 0 or len(target) == 0:
        return []

    tree = cKDTree(source_heavy.positions)
    keys: list[ResidueKey] = []
    for residue in target.residues:
        heavy = residue.atoms.select_atoms("not name H*")
        if len(heavy) == 0:
            continue
        if any(tree.query_ball_point(heavy.positions, cutoff_A)):
            keys.append(
                ResidueKey(
                    chain=target_chain,
                    resid=residue_number_with_icode(residue),
                    resname=str(residue.resname).strip() or None,
                )
            )
    return keys


def parse_biolip_record(line: str) -> BioLiPRecord:
    """Parse the columns needed from one BioLiP/BioLiP2 annotation row.

    The BioLiP flat-file format is tab-separated. For MolInspect static spatial
    validation, the essential columns are PDB ID, receptor chain, ligand ID,
    ligand chain, ligand serial number, and binding-site residue tokens.
    """

    columns = line.rstrip("\n").split("\t")
    if len(columns) < 9:
        raise ValueError("BioLiP row has fewer than 9 tab-separated columns")
    receptor_chain = columns[1].strip()
    binding_site = columns[7].strip()
    ligand_auth_resid = _clean_optional(columns[19]) if len(columns) > 19 else None
    return BioLiPRecord(
        pdb_id=columns[0].strip().lower(),
        receptor_chain=receptor_chain,
        resolution_A=_as_float(columns[2]),
        binding_site_code=_clean_optional(columns[3]),
        ligand_id=columns[4].strip(),
        ligand_chain=columns[5].strip(),
        ligand_serial=columns[6].strip(),
        ligand_auth_resid=ligand_auth_resid,
        binding_site_residues=_parse_biolip_binding_site(binding_site, receptor_chain),
    )


def parse_p2rank_dataset_items(text: str) -> list[str]:
    """Parse P2Rank `.ds` text into dataset item paths."""

    items: list[str] = []
    for line in text.splitlines():
        clean = line.split("#", 1)[0].strip()
        if not clean:
            continue
        items.extend(part for part in clean.split() if part)
    return items


def _parse_biolip_binding_site(binding_site: str, chain: str) -> list[ResidueKey]:
    keys: list[ResidueKey] = []
    for token in binding_site.split():
        resid = _biolip_resid_from_token(token)
        if resid:
            keys.append(ResidueKey(chain=chain, resid=resid))
    return keys


def _biolip_resid_from_token(token: str) -> str | None:
    token = token.strip()
    if len(token) < 2:
        return None
    # BioLiP residue tokens encode the amino-acid one-letter code followed by the
    # residue number/insertion token, for example H87 or D189A.
    if token[0].isalpha() and any(char.isdigit() for char in token[1:]):
        return token[1:]
    return token if any(char.isdigit() for char in token) else None


def _safe_ratio(numerator: float, denominator: float) -> float:
    if denominator == 0:
        return 1.0 if numerator == 0 else 0.0
    return round(numerator / denominator, 6)


def _as_float(value: str) -> float | None:
    try:
        parsed = float(value)
    except ValueError:
        return None
    return parsed if parsed >= 0 else None


def _clean_optional(value: str) -> str | None:
    clean = value.strip()
    return clean or None


def _protein_chain_atomgroup(universe: Any, chain: str) -> Any:
    try:
        atoms = universe.select_atoms(f"protein and chainID {chain}")
    except Exception:
        atoms = universe.atoms[:0]
    if len(atoms) > 0:
        return atoms
    try:
        return universe.select_atoms(f"protein and segid {chain}")
    except Exception:
        return universe.atoms[:0]
