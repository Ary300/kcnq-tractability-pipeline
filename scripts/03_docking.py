#!/usr/bin/env python3
from __future__ import annotations

import sys
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from kcnq_pipeline.config import EXPORT_DIR
from kcnq_pipeline.docking import prepare_ligands, prepare_receptors, run_full_matrix_docking, run_top10_docking
from kcnq_pipeline.structures import fetch_structures, parse_structure_coords


def main() -> None:
    top10 = pd.read_csv(EXPORT_DIR / "top10_scored.csv")
    structure_files = fetch_structures()
    coords, centroids = parse_structure_coords(structure_files)
    ligands = prepare_ligands(EXPORT_DIR / "ligands")
    receptors, boxes = prepare_receptors(structure_files, coords, centroids, EXPORT_DIR / "receptors")
    top10_docking = run_top10_docking(top10, receptors, boxes, ligands)
    matrix = run_full_matrix_docking(top10, receptors, boxes, ligands)
    top10_docking.to_csv(EXPORT_DIR / "top10_docking.csv", index=False)
    matrix.to_csv(EXPORT_DIR / "docking_matrix.csv", index=False)
    print("Top10 docked:", len(top10_docking))
    print("Matrix rows:", len(matrix))


if __name__ == "__main__":
    main()
