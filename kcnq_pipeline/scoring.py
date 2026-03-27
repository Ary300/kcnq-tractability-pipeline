from __future__ import annotations

import io
import json
import math
from pathlib import Path

import numpy as np
import pandas as pd
import requests

from .config import CACHE_DIR, CADD_URL, CLINVAR_TIER_MAP, EXCLUDE_SIG, REGION_DRUG_MAP
from .utils import load_json, rarity_score, read_gz_tsv_subset, session, wait_between_requests


def add_clinvar_and_rarity_scores(df: pd.DataFrame) -> pd.DataFrame:
    out = df.copy()
    out["clinvar_tier_score"] = out["clinical_significance"].map(lambda x: CLINVAR_TIER_MAP.get(x, 0.2))
    out["rarity_score"] = out["gnomad_af"].apply(rarity_score) if "gnomad_af" in out.columns else 1.0
    out["region_drug_score"] = out["region_label"].map(REGION_DRUG_MAP).fillna(0.2)
    return out


def fetch_cadd_scores(df: pd.DataFrame) -> pd.DataFrame:
    out = df.copy()
    cache_path = CACHE_DIR / "cadd_cache.json"
    cache = load_json(cache_path, default={}) or {}
    variant_df = out[["chromosome", "position", "ref", "alt"]].dropna().copy()
    variant_df = variant_df[
        variant_df["ref"].astype(str).str.len().eq(1) & variant_df["alt"].astype(str).str.len().eq(1)
    ].copy()
    variant_df["position"] = variant_df["position"].astype(int)
    variant_keys = sorted(
        {
            f"{row.chromosome}:{int(row.position)}:{row.ref}:{row.alt}"
            for row in variant_df.itertuples(index=False)
        }
    )
    missing = [key for key in variant_keys if key not in cache]
    s = session()
    if missing:
        missing_df = pd.DataFrame(
            [key.split(":") for key in missing],
            columns=["chromosome", "position", "ref", "alt"],
        )
        missing_df["position"] = missing_df["position"].astype(int)
        windows: list[tuple[str, int, int]] = []
        for chrom, grp in missing_df[["chromosome", "position"]].drop_duplicates().sort_values(
            ["chromosome", "position"]
        ).groupby("chromosome"):
            positions = grp["position"].tolist()
            start = positions[0]
            end = positions[0]
            for pos in positions[1:]:
                if pos - start <= 99:
                    end = pos
                else:
                    windows.append((chrom, start, end))
                    start = pos
                    end = pos
            windows.append((chrom, start, end))
    else:
        windows = []
    missing_set = set(missing)
    for idx, (chrom, start, end) in enumerate(windows, start=1):
        url = f"{CADD_URL}/{chrom}:{start}-{end}"
        try:
            resp = s.get(url, timeout=30)
            if resp.status_code == 200:
                data = resp.json()
                rows = data[1:] if isinstance(data, list) and data and isinstance(data[0], list) else []
                for row in rows:
                    key = f"{row[0]}:{int(row[1])}:{row[2]}:{row[3]}"
                    if key in missing_set:
                        cache[key] = {"cadd_phred": row[5], "cadd_raw": row[4]}
        except Exception:
            pass
        wait_between_requests(0.05)
        if idx % 25 == 0 or idx == len(windows):
            print(f"CADD windows: {idx}/{len(windows)} queried")
    for key in missing:
        cache.setdefault(key, {"cadd_phred": np.nan, "cadd_raw": np.nan})
    cache_path.write_text(json.dumps(cache, indent=2))
    out["variant_key"] = out.apply(
        lambda r: f"{r['chromosome']}:{int(r['position'])}:{r['ref']}:{r['alt']}"
        if pd.notna(r["position"]) and pd.notna(r["ref"]) and pd.notna(r["alt"])
        else pd.NA,
        axis=1,
    )
    out["cadd_phred"] = out["variant_key"].map(lambda k: cache.get(k, {}).get("cadd_phred") if pd.notna(k) else np.nan)
    out["cadd_norm"] = out["cadd_phred"].apply(lambda x: min(max(float(x) / 40.0, 0.0), 1.0) if pd.notna(x) else np.nan)
    return out


def fetch_alphamissense_scores(df: pd.DataFrame, alphamissense_path: Path) -> pd.DataFrame:
    out = df.copy()
    out["chromosome"] = out["chromosome"].astype("string")
    out["position"] = out["position"].astype("Int64")
    out["ref"] = out["ref"].astype("string")
    out["alt"] = out["alt"].astype("string")
    keys = {
        (str(chrom), int(pos), str(ref), str(alt))
        for chrom, pos, ref, alt in out[["chromosome", "position", "ref", "alt"]].dropna().itertuples(index=False)
    }
    print(f"AlphaMissense: scanning {len(keys)} keys from {alphamissense_path.name}")
    subset = read_gz_tsv_subset(alphamissense_path, keys)
    subset["chromosome"] = subset["#CHROM"].str.replace("chr", "", regex=False)
    subset["position"] = subset["POS"].astype(int)
    subset = subset.rename(
        columns={"REF": "ref", "ALT": "alt", "am_pathogenicity": "alphamissense_score", "am_class": "alphamissense_class"}
    )
    subset["chromosome"] = subset["chromosome"].astype("string")
    subset["position"] = subset["position"].astype("Int64")
    subset["ref"] = subset["ref"].astype("string")
    subset["alt"] = subset["alt"].astype("string")
    subset = subset[["chromosome", "position", "ref", "alt", "alphamissense_score", "alphamissense_class"]].drop_duplicates()
    print(f"AlphaMissense: matched {len(subset)} variants")
    out = out.merge(subset, on=["chromosome", "position", "ref", "alt"], how="left")
    return out


def apply_composite_pathogenicity(df: pd.DataFrame) -> pd.DataFrame:
    out = df.copy()
    full_signal = out["cadd_norm"].notna() & out["alphamissense_score"].notna()
    out["path_score_imputed"] = ~full_signal
    out["path_score"] = np.where(
        full_signal,
        0.40 * out["clinvar_tier_score"] + 0.30 * out["alphamissense_score"] + 0.30 * out["cadd_norm"],
        out["clinvar_tier_score"],
    )
    return out


def pocket_proximity_score(dist: float | None) -> float:
    if dist is None or pd.isna(dist):
        return 0.3
    return float(math.exp(-float(dist) / 15.0))


def apply_rescue_scores(df: pd.DataFrame) -> pd.DataFrame:
    out = df.copy()
    out["pocket_proximity_score"] = out["pocket_dist_A"].apply(pocket_proximity_score)
    out["rescue_score"] = (
        0.30 * out["path_score"]
        + 0.35 * out["pocket_proximity_score"]
        + 0.25 * out["region_drug_score"]
        + 0.10 * out["rarity_score"]
    )
    return out


def _tractability_assignment(row: pd.Series) -> tuple[float, str, str]:
    if row.get("gene") == "KCNQ2":
        rescue_class = row.get("bhatt_wt_cotransfection_rescue")
        if pd.notna(rescue_class) and rescue_class == "-D N-":
            return 0.50, "dominant_negative", "bhatt_direct"
        if pd.notna(rescue_class) and rescue_class == "-f R-":
            return 1.00, "full_rescue", "bhatt_direct"
        if pd.notna(rescue_class) and rescue_class == "-p R-":
            return 1.00, "partial_rescue", "bhatt_direct"
        if pd.notna(rescue_class) and rescue_class == "none":
            return 0.75, "no_rescue_annotated", "bhatt_direct"
        if (
            row.get("region_label") == "Voltage-Sensing Domain"
            and pd.notna(row.get("bhatt_homo_current_rel_wt"))
            and float(row["bhatt_homo_current_rel_wt"]) > 0.30
        ):
            return 0.95, "vsd_retained_current", "bhatt_proxy"
        if (
            row.get("region_label") in {"Pore Domain (S5-S6)", "Pore/Selectivity Filter"}
            and pd.notna(row.get("bhatt_homo_current_rel_wt"))
            and float(row["bhatt_homo_current_rel_wt"]) <= 0.10
        ):
            return 0.75, "severe_pore_lof", "bhatt_proxy"
    if row.get("region_label") == "Voltage-Sensing Domain":
        return 0.95, "vsd_gating_like", "heuristic"
    if (
        row.get("region_label") in {"Pore Domain (S5-S6)", "Pore/Selectivity Filter"}
        and pd.notna(row.get("cadd_phred"))
        and float(row["cadd_phred"]) >= 30.0
    ):
        return 0.75, "severe_pore_lof_like", "heuristic"
    if row.get("region_label") in {"CaM/Interface Region", "Trafficking/Assembly"}:
        return 0.80, "trafficking_interface_like", "heuristic"
    return 0.85, "default_unknown", "heuristic"


def apply_tractability_modifier(df: pd.DataFrame, bhatt_df: pd.DataFrame | None = None) -> pd.DataFrame:
    out = df.copy()
    if bhatt_df is not None and not bhatt_df.empty:
        bhatt_cols = [
            "protein_change",
            "bhatt_homo_current_rel_wt",
            "bhatt_hetero_current_rel_wt",
            "bhatt_homo_functional_class",
            "bhatt_hetero_functional_class",
            "bhatt_wt_cotransfection_rescue",
        ]
        bhatt_merge = bhatt_df[bhatt_cols].drop_duplicates("protein_change").copy()
        bhatt_merge["bhatt_benchmark_gene"] = "KCNQ2"
        out = out.merge(
            bhatt_merge,
            left_on=["gene", "protein_change"],
            right_on=["bhatt_benchmark_gene", "protein_change"],
            how="left",
        )
        out = out.drop(columns=["bhatt_benchmark_gene"])
    out["structural_opportunity_score"] = out["rescue_score"]
    tractability = out.apply(_tractability_assignment, axis=1, result_type="expand")
    tractability.columns = ["tractability_modifier", "tractability_class", "tractability_basis"]
    out = pd.concat([out, tractability], axis=1)
    out["rescue_priority_score"] = out["structural_opportunity_score"] * out["tractability_modifier"]
    return out


def clinical_subset(df: pd.DataFrame) -> pd.DataFrame:
    return df[(df["source"] == "ClinVar") & (~df["clinical_significance"].isin(EXCLUDE_SIG))].copy()


def select_top10(df_clinical: pd.DataFrame) -> pd.DataFrame:
    score_col = "rescue_priority_score" if "rescue_priority_score" in df_clinical.columns else "rescue_score"
    ranked = df_clinical.sort_values(score_col, ascending=False).reset_index(drop=True)
    selected = []
    counts: dict[str, int] = {}
    for _, row in ranked.iterrows():
        gene = row["gene"]
        if counts.get(gene, 0) >= 3:
            continue
        selected.append(row)
        counts[gene] = counts.get(gene, 0) + 1
        if len(selected) == 10:
            break
    out = pd.DataFrame(selected).reset_index(drop=True)
    out["priority_rank"] = out.index + 1
    return out


def assign_candidate_drug(df: pd.DataFrame, drugs_df: pd.DataFrame) -> pd.DataFrame:
    out = df.copy()
    names = []
    phases = []
    for _, row in out.iterrows():
        gene_drugs = drugs_df[drugs_df["gene"] == row["gene"]]
        if gene_drugs.empty:
            names.append("None")
            phases.append(0)
            continue
        matched = gene_drugs[gene_drugs["binding_region"] == row["region_label"]]
        source = matched if not matched.empty else gene_drugs
        best = source.sort_values("max_phase", ascending=False).iloc[0]
        names.append(best["drug_name"])
        phases.append(best["max_phase"])
    out["candidate_drug"] = names
    out["drug_phase"] = phases
    return out
