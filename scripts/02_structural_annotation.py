from __future__ import annotations

import sys
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from kcnq_pipeline.config import EXPORT_DIR
from kcnq_pipeline.structures import annotate_variants, fetch_structures, parse_structure_coords


def main() -> None:
    merged = pd.read_csv(EXPORT_DIR / "all_variants_merged.csv")
    structure_files = fetch_structures()
    coords, centroids = parse_structure_coords(structure_files)
    annotated = annotate_variants(merged, coords, centroids)
    annotated.to_csv(EXPORT_DIR / "all_variants_structural.csv", index=False)
    print("Annotated:", len(annotated))
    print("With coords:", int(annotated["structural_coords_available"].sum()))


if __name__ == "__main__":
    main()
