"""Fetch small real PDB examples from RCSB.

The examples are intentionally modest and canonical. They are used for local validation
and visual QA, not as the primary offline unit-test fixtures.
"""

from __future__ import annotations

from pathlib import Path
from urllib.request import urlopen

EXAMPLES = {
    "1CRN": "Tiny single-chain crambin protein.",
    "1UBQ": "Canonical compact ubiquitin protein.",
    "4HHB": "Classic hemoglobin tetramer with heme ligands.",
}

BASE_URL = "https://files.rcsb.org/download/{pdb_id}.pdb"
DATA_DIR = Path(__file__).resolve().parents[1] / "examples" / "data"


def main() -> None:
    DATA_DIR.mkdir(parents=True, exist_ok=True)
    for pdb_id, description in EXAMPLES.items():
        url = BASE_URL.format(pdb_id=pdb_id)
        target = DATA_DIR / f"{pdb_id.lower()}.pdb"
        print(f"Fetching {pdb_id}: {description}")
        with urlopen(url, timeout=30) as response:
            content = response.read()
        if not content.startswith((b"HEADER", b"TITLE", b"COMPND")):
            raise RuntimeError(f"Unexpected response for {pdb_id} from {url}")
        target.write_bytes(content)
        print(f"  wrote {target} ({len(content):,} bytes)")


if __name__ == "__main__":
    main()
