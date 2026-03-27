from __future__ import annotations

import sys
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from kcnq_pipeline.config import ALPHAMISSENSE_URL, BENCHMARK_DIR, EXPORT_DIR, RAW_DIR, KNOWN_KCNQ_DRUGS
from kcnq_pipeline.benchmarking import parse_bhatt_kcnq2
from kcnq_pipeline.scoring import (
    add_clinvar_and_rarity_scores,
    apply_composite_pathogenicity,
    apply_rescue_scores,
    apply_tractability_modifier,
    assign_candidate_drug,
    clinical_subset,
    fetch_alphamissense_scores,
    fetch_cadd_scores,
    select_top10,
)
from kcnq_pipeline.utils import download


def main() -> None:
    annotated = pd.read_csv(EXPORT_DIR / "all_variants_structural.csv")
    df = add_clinvar_and_rarity_scores(annotated)
    df = fetch_cadd_scores(df)
    download(ALPHAMISSENSE_URL, RAW_DIR / "AlphaMissense_hg38.tsv.gz")
    am_path = RAW_DIR / "AlphaMissense_hg38.tsv.gz"
    df = fetch_alphamissense_scores(df, am_path)
    df = apply_composite_pathogenicity(df)
    df = apply_rescue_scores(df)
    bhatt = parse_bhatt_kcnq2() if (BENCHMARK_DIR / "bhatt_kcnq2_supp_table4.xlsx").exists() else pd.DataFrame()
    df = apply_tractability_modifier(df, bhatt_df=bhatt)
    clinical = clinical_subset(df)
    top10 = select_top10(clinical)
    drugs_df = pd.DataFrame(KNOWN_KCNQ_DRUGS)
    top10 = assign_candidate_drug(top10, drugs_df)
    df.to_csv(EXPORT_DIR / "all_variants_scored.csv", index=False)
    clinical.to_csv(EXPORT_DIR / "clinical_variants_scored.csv", index=False)
    top10.to_csv(EXPORT_DIR / "top10_scored.csv", index=False)
    print("Merged scored:", len(df))
    print("Clinical ranked:", len(clinical))
    print("Top10:", len(top10))
    print("Full CADD+AM coverage:", int((df['cadd_norm'].notna() & df['alphamissense_score'].notna()).sum()))
    print("Imputed path score:", int(df["path_score_imputed"].sum()))
    print("Bhatt direct tractability rows:", int((df["tractability_basis"] == "bhatt_direct").sum()))


if __name__ == "__main__":
    main()
