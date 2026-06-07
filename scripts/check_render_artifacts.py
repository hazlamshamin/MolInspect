"""Check generated render PNGs are present and nonblank."""

from __future__ import annotations

from pathlib import Path

import numpy as np
from PIL import Image

ROOT = Path(__file__).resolve().parents[1]
RENDER_DIR = ROOT / "examples" / "renders"
EXPECTED = ["1crn_overview.png", "1ubq_overview.png", "4hhb_overview.png"]
EXPECTED += ["1crn_pymol.png", "1ubq_pymol.png", "4hhb_pymol.png"]


def main() -> None:
    for filename in EXPECTED:
        path = RENDER_DIR / filename
        if not path.exists():
            raise FileNotFoundError(path)

        image = Image.open(path).convert("RGB")
        array = np.asarray(image)
        nonwhite_fraction = np.mean(np.any(array < 245, axis=2))
        if image.size[0] < 640 or image.size[1] < 480:
            raise RuntimeError(f"{path} is unexpectedly small: {image.size}")
        if nonwhite_fraction < 0.02:
            raise RuntimeError(f"{path} appears blank: nonwhite_fraction={nonwhite_fraction:.4f}")
        print(f"{path.name}: size={image.size}, nonwhite_fraction={nonwhite_fraction:.4f}")


if __name__ == "__main__":
    main()
