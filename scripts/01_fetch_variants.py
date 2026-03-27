from __future__ import annotations

import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from kcnq_pipeline.config import EXPORT_DIR
from kcnq_pipeline.fetch import build_variant_tables


def main() -> None:
    clinvar, gnomad, merged = build_variant_tables()
    clinvar.to_csv(EXPORT_DIR / "clinvar_uncapped.csv", index=False)
    gnomad.to_csv(EXPORT_DIR / "gnomad_uncapped.csv", index=False)
    merged.to_csv(EXPORT_DIR / "all_variants_merged.csv", index=False)
    print("ClinVar:", len(clinvar))
    print("gnomAD:", len(gnomad))
    print("Merged dedup:", len(merged))


if __name__ == "__main__":
    main()
