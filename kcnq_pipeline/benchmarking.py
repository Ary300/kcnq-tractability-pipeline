from __future__ import annotations

import itertools
from pathlib import Path
import subprocess

import numpy as np
import pandas as pd
from scipy import stats
import statsmodels.formula.api as smf
from Bio import AlignIO, SeqIO
from sklearn.metrics import average_precision_score, roc_auc_score
from sklearn.preprocessing import StandardScaler

from .config import BENCHMARK_DIR, CACHE_DIR, UNIPROT_IDS
from .utils import download, session


BASE_COMPONENT_WEIGHTS = {
    "path": 0.30,
    "pocket": 0.35,
    "region": 0.25,
    "rarity": 0.10,
}

ORTHOLOG_SPECIES = {
    "human": 9606,
    "mouse": 10090,
    "rat": 10116,
    "zebrafish": 7955,
}

AA_GROUPS = {
    "A": "hydrophobic",
    "V": "hydrophobic",
    "L": "hydrophobic",
    "I": "hydrophobic",
    "M": "hydrophobic",
    "F": "aromatic",
    "W": "aromatic",
    "Y": "aromatic",
    "S": "polar",
    "T": "polar",
    "N": "polar",
    "Q": "polar",
    "K": "basic",
    "R": "basic",
    "H": "basic",
    "D": "acidic",
    "E": "acidic",
    "C": "special",
    "G": "special",
    "P": "special",
}


def bootstrap_median_difference(a: pd.Series, b: pd.Series, n_boot: int = 2000) -> tuple[float, float, float]:
    a = a.dropna().to_numpy()
    b = b.dropna().to_numpy()
    rng = np.random.default_rng(42)
    diffs = []
    for _ in range(n_boot):
        sa = rng.choice(a, size=len(a), replace=True)
        sb = rng.choice(b, size=len(b), replace=True)
        diffs.append(np.median(sa) - np.median(sb))
    lo, hi = np.percentile(diffs, [2.5, 97.5])
    return float(np.median(a) - np.median(b)), float(lo), float(hi)


def cohens_d(a: pd.Series, b: pd.Series) -> float:
    a = a.dropna().to_numpy()
    b = b.dropna().to_numpy()
    pooled = np.sqrt(((len(a) - 1) * np.var(a, ddof=1) + (len(b) - 1) * np.var(b, ddof=1)) / (len(a) + len(b) - 2))
    return float((np.mean(a) - np.mean(b)) / pooled) if pooled else float("nan")


def section32_stats(df_clinical: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    tmp = df_clinical.copy()
    tmp["sig_group"] = np.where(
        tmp["clinical_significance"].isin(["Pathogenic", "Likely pathogenic", "Pathogenic/Likely pathogenic"]),
        "P/LP",
        np.where(tmp["clinical_significance"] == "Uncertain significance", "VUS", pd.NA),
    )
    tmp = tmp[tmp["sig_group"].isin(["P/LP", "VUS"])].copy()
    tmp["pathogenic_binary"] = (tmp["sig_group"] == "P/LP").astype(int)
    plp = tmp[tmp["pathogenic_binary"] == 1]["pocket_dist_A"]
    vus = tmp[tmp["pathogenic_binary"] == 0]["pocket_dist_A"]
    u_stat, p_value = stats.mannwhitneyu(plp.dropna(), vus.dropna(), alternative="two-sided")
    med_diff, ci_lo, ci_hi = bootstrap_median_difference(plp, vus, n_boot=10_000)
    effect = cohens_d(plp, vus)
    summary = pd.DataFrame(
        [
            {
                "comparison": "P/LP vs VUS pocket distance",
                "n_plp": int(plp.dropna().shape[0]),
                "n_vus": int(vus.dropna().shape[0]),
                "mannwhitney_u": u_stat,
                "p_value": p_value,
                "cohens_d": effect,
                "median_difference": med_diff,
                "median_diff_ci_lo": ci_lo,
                "median_diff_ci_hi": ci_hi,
            }
        ]
    )
    reg = tmp.dropna(subset=["pocket_dist_A", "gene"]).copy()
    model = smf.logit("pathogenic_binary ~ pocket_dist_A + C(gene)", data=reg).fit(disp=False, maxiter=1000)
    coeffs = model.summary2().tables[1].reset_index().rename(columns={"index": "term"})
    coeffs["odds_ratio"] = np.exp(coeffs["Coef."])
    coeffs["or_ci_lo"] = np.exp(coeffs["Coef."] - 1.96 * coeffs["Std.Err."])
    coeffs["or_ci_hi"] = np.exp(coeffs["Coef."] + 1.96 * coeffs["Std.Err."])
    return summary, coeffs


def _structural_score_from_weights(df: pd.DataFrame, weights: dict[str, float]) -> pd.Series:
    return (
        weights["path"] * df["path_score"]
        + weights["pocket"] * df["pocket_proximity_score"]
        + weights["region"] * df["region_drug_score"]
        + weights["rarity"] * df["rarity_score"]
    )


def _select_diverse_top10(df: pd.DataFrame, score_col: str) -> pd.DataFrame:
    ranked = df.sort_values(score_col, ascending=False).reset_index(drop=True)
    chosen = []
    counts: dict[str, int] = {}
    for _, row in ranked.iterrows():
        gene = row["gene"]
        if counts.get(gene, 0) >= 3:
            continue
        chosen.append(row)
        counts[gene] = counts.get(gene, 0) + 1
        if len(chosen) == 10:
            break
    out = pd.DataFrame(chosen).reset_index(drop=True)
    out["trial_rank"] = out.index + 1
    return out


def weight_sensitivity(df_clinical: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    labels = ["path", "pocket", "region", "rarity"]
    deltas = [-0.20, -0.10, 0.0, 0.10, 0.20]
    base_top = _select_diverse_top10(df_clinical, "rescue_priority_score" if "rescue_priority_score" in df_clinical.columns else "rescue_score")
    base_labels = base_top.apply(lambda r: f"{r['gene']} {r['protein_change']}", axis=1).tolist()
    grid_rows = []
    rank_rows = []
    for combo_id, factors in enumerate(itertools.product(deltas, repeat=4), start=1):
        trial = {
            label: BASE_COMPONENT_WEIGHTS[label] * (1.0 + delta)
            for label, delta in zip(labels, factors)
        }
        total = sum(trial.values())
        trial = {k: v / total for k, v in trial.items()}
        tmp = df_clinical.copy()
        tmp["structural_score_trial"] = _structural_score_from_weights(tmp, trial)
        if "tractability_modifier" in tmp.columns:
            tmp["priority_score_trial"] = tmp["structural_score_trial"] * tmp["tractability_modifier"]
            score_col = "priority_score_trial"
        else:
            score_col = "structural_score_trial"
        top = _select_diverse_top10(tmp, score_col)
        labels_top = top.apply(lambda r: f"{r['gene']} {r['protein_change']}", axis=1).tolist()
        grid_rows.append(
            {
                "combo_id": combo_id,
                "path_delta": factors[0],
                "pocket_delta": factors[1],
                "region_delta": factors[2],
                "rarity_delta": factors[3],
                "path_weight": trial["path"],
                "pocket_weight": trial["pocket"],
                "region_weight": trial["region"],
                "rarity_weight": trial["rarity"],
                "top10_overlap_with_base": len(set(base_labels) & set(labels_top)),
                "membership_changed": len(set(base_labels) ^ set(labels_top)),
                "top10_labels": " | ".join(labels_top),
            }
        )
        for _, row in top.iterrows():
            rank_rows.append(
                {
                    "combo_id": combo_id,
                    "candidate": f"{row['gene']} {row['protein_change']}",
                    "trial_rank": int(row["trial_rank"]),
                    "path_delta": factors[0],
                    "pocket_delta": factors[1],
                    "region_delta": factors[2],
                    "rarity_delta": factors[3],
                }
            )
    grid_df = pd.DataFrame(grid_rows)
    rank_df = pd.DataFrame(rank_rows)
    total_combos = len(grid_df)

    freq_rows = []
    observed_candidates = sorted(rank_df["candidate"].unique())
    for candidate in observed_candidates:
        subset = rank_df[rank_df["candidate"] == candidate]
        freq_rows.append(
            {
                "candidate": candidate,
                "is_base_top10": candidate in base_labels,
                "inclusion_count": int(len(subset)),
                "inclusion_proportion": float(len(subset) / total_combos),
                "present_in_ge_90pct": float(len(subset) / total_combos) >= 0.90,
                "median_rank_when_present": float(subset["trial_rank"].median()),
                "mean_rank_when_present": float(subset["trial_rank"].mean()),
                "min_rank_when_present": int(subset["trial_rank"].min()),
                "max_rank_when_present": int(subset["trial_rank"].max()),
            }
        )
    freq_df = pd.DataFrame(freq_rows).sort_values(
        ["present_in_ge_90pct", "inclusion_proportion", "median_rank_when_present"],
        ascending=[False, False, True],
    )

    effect_rows = []
    for candidate in observed_candidates:
        cand = rank_df[rank_df["candidate"] == candidate].copy()
        for weight_name in labels:
            delta_col = f"{weight_name}_delta"
            grouped = cand.groupby(delta_col)["trial_rank"].agg(["count", "median"]).reset_index()
            presence = grouped.set_index(delta_col)["count"].to_dict()
            medians = grouped.set_index(delta_col)["median"].to_dict()
            effect_rows.append(
                {
                    "candidate": candidate,
                    "weight_name": weight_name,
                    "presence_rate_span": (max(presence.values()) - min(presence.values())) / total_combos if presence else 0.0,
                    "median_rank_span": max(medians.values()) - min(medians.values()) if medians else np.nan,
                    "median_rank_at_neg20": medians.get(-0.20, np.nan),
                    "median_rank_at_neg10": medians.get(-0.10, np.nan),
                    "median_rank_at_base": medians.get(0.0, np.nan),
                    "median_rank_at_pos10": medians.get(0.10, np.nan),
                    "median_rank_at_pos20": medians.get(0.20, np.nan),
                }
            )
    effect_df = pd.DataFrame(effect_rows).sort_values(
        ["presence_rate_span", "median_rank_span"], ascending=[False, False]
    )
    return grid_df, freq_df, effect_df


def leave_one_component_out_auc(df_clinical: pd.DataFrame, bhatt_df: pd.DataFrame) -> pd.DataFrame:
    kcnq2 = df_clinical[df_clinical["gene"] == "KCNQ2"].copy()
    overlap = kcnq2.merge(
        annotate_bhatt_mechanism(bhatt_df.copy()),
        on="protein_change",
        how="inner",
    )
    rows = []
    models = {"full_model": BASE_COMPONENT_WEIGHTS.copy()}
    for dropped in ["path", "pocket", "region", "rarity"]:
        trial = {k: v for k, v in BASE_COMPONENT_WEIGHTS.items() if k != dropped}
        total = sum(trial.values())
        trial = {k: v / total for k, v in trial.items()}
        trial[dropped] = 0.0
        models[f"drop_{dropped}"] = trial
    for model_name, weights in models.items():
        scored = overlap.copy()
        scored["structural_score_trial"] = _structural_score_from_weights(scored, weights)
        scored["priority_score_trial"] = scored["structural_score_trial"] * scored["tractability_modifier"]
        rows.append(
            {
                "model_name": model_name,
                "path_weight": weights["path"],
                "pocket_weight": weights["pocket"],
                "region_weight": weights["region"],
                "rarity_weight": weights["rarity"],
                "auc_structural_responsive_current50": _score_endpoint_auc(scored, "structural_score_trial", "bhatt_retig_responsive_current50"),
                "auc_structural_responsive_any": _score_endpoint_auc(scored, "structural_score_trial", "bhatt_retig_responsive_any"),
                "auc_structural_reaches_wt_current": _score_endpoint_auc(scored, "structural_score_trial", "bhatt_retig_reaches_wt_current"),
                "auc_priority_responsive_current50": _score_endpoint_auc(scored, "priority_score_trial", "bhatt_retig_responsive_current50"),
                "auc_priority_responsive_any": _score_endpoint_auc(scored, "priority_score_trial", "bhatt_retig_responsive_any"),
                "auc_priority_reaches_wt_current": _score_endpoint_auc(scored, "priority_score_trial", "bhatt_retig_reaches_wt_current"),
            }
        )
    return pd.DataFrame(rows).sort_values("auc_structural_responsive_any", ascending=False).reset_index(drop=True)


def build_kcnq4_literature_table(df_clinical: pd.DataFrame, df_top10: pd.DataFrame) -> pd.DataFrame:
    curated = {
        "L281S": {
            "published_homomeric_current_pct_wt": "Near-zero / severely reduced",
            "published_opener_response_homomeric": "No meaningful rescue reported for homomeric channels",
            "published_opener_response_heteromeric": "ZnP + retigabine restored WT/L281S current toward near-WT",
            "surface_expression": "Reduced surface expression",
            "mechanism_class": "Dual mechanism: trafficking defect plus conductance loss; heteromeric rescue relevant",
            "clinical_phenotype_literature": "DFNA2 / nonsyndromic hearing loss",
            "source": "Leitner 2012; Gao 2013",
            "literature_status": "Direct functional support",
        },
        "T280K": {
            "published_homomeric_current_pct_wt": "No direct report found",
            "published_opener_response_homomeric": "No direct report found",
            "published_opener_response_heteromeric": "No direct report found",
            "surface_expression": "No direct report found",
            "mechanism_class": "Unknown; current pipeline prediction only",
            "clinical_phenotype_literature": "No exact PubMed title/abstract hit for KCNQ4 T280K / Thr280Lys as of 2026-03-22",
            "source": "PubMed exact-variant search on 2026-03-22 returned no items",
            "literature_status": "Novel / uncataloged candidate",
        },
        "G285C": {
            "published_homomeric_current_pct_wt": "Severely reduced",
            "published_opener_response_homomeric": "No variant-specific opener rescue report found",
            "published_opener_response_heteromeric": "No variant-specific opener rescue report found",
            "surface_expression": "Reduced surface expression",
            "mechanism_class": "Trafficking and conductance defect",
            "clinical_phenotype_literature": "DFNA2 / nonsyndromic hearing loss",
            "source": "Gao 2013",
            "literature_status": "Direct functional support",
        },
        "G285S": {
            "published_homomeric_current_pct_wt": "Near-zero / absent",
            "published_opener_response_homomeric": "No meaningful rescue reported for homomeric channels",
            "published_opener_response_heteromeric": "WT/G285S current rescued strongly by ZnP + retigabine",
            "surface_expression": "Reduced surface expression",
            "mechanism_class": "Severe pore LOF in homomers; heteromeric dominant-negative rescue context",
            "clinical_phenotype_literature": "DFNA2 / nonsyndromic hearing loss",
            "source": "Leitner 2012; Gao 2013; Lee 2021",
            "literature_status": "Direct functional support",
        },
    }
    current_top = set(df_top10.loc[df_top10["gene"] == "KCNQ4", "protein_change"].tolist())
    source_df = df_clinical[(df_clinical["gene"] == "KCNQ4") & (df_clinical["protein_change"].isin(curated))].copy()
    rows = []
    for variant, notes in curated.items():
        row_df = source_df[source_df["protein_change"] == variant]
        row = row_df.iloc[0].to_dict() if not row_df.empty else {}
        rows.append(
            {
                "variant": f"KCNQ4 {variant}",
                "in_current_top10": variant in current_top,
                "structural_opportunity_score": row.get("structural_opportunity_score", row.get("rescue_score", np.nan)),
                "tractability_modifier": row.get("tractability_modifier", np.nan),
                "rescue_priority_score": row.get("rescue_priority_score", np.nan),
                "clinvar_significance": row.get("clinical_significance", pd.NA),
                "clinvar_trait": row.get("trait", pd.NA),
                "published_homomeric_current_pct_wt": notes["published_homomeric_current_pct_wt"],
                "published_opener_response_homomeric": notes["published_opener_response_homomeric"],
                "published_opener_response_heteromeric": notes["published_opener_response_heteromeric"],
                "surface_expression": notes["surface_expression"],
                "mechanism_class": notes["mechanism_class"],
                "clinical_phenotype_literature": notes["clinical_phenotype_literature"],
                "literature_status": notes["literature_status"],
                "source": notes["source"],
            }
        )
    return pd.DataFrame(rows)


def parse_bhatt_kcnq2() -> pd.DataFrame:
    path = BENCHMARK_DIR / "bhatt_kcnq2_supp_table4.xlsx"
    rel = pd.read_excel(path, sheet_name="Relative to WT channel", header=None)
    rel = rel.iloc[16:].copy()
    rel = rel[rel[3].notna()].copy()
    rel["protein_change"] = rel[3].astype(str)
    rel["bhatt_current_rel_wt"] = pd.to_numeric(rel[8], errors="coerce")
    rel["bhatt_delta_vhalf"] = pd.to_numeric(rel[14], errors="coerce")
    rel["bhatt_k_ratio"] = pd.to_numeric(rel[18], errors="coerce")
    rel["bhatt_tau_ratio"] = pd.to_numeric(rel[23], errors="coerce")
    rel["bhatt_phenotype"] = rel[5].astype("string")

    func = pd.read_excel(path, sheet_name="Functional phenotype & rescue", header=None)
    func = func.iloc[24:].copy()
    func = func[func[3].notna()].copy()
    func["protein_change"] = func[3].astype(str)
    func["bhatt_homo_current_rel_wt"] = pd.to_numeric(func[11], errors="coerce")
    func["bhatt_homo_functional_class"] = func[13].astype("string")
    func["bhatt_hetero_current_rel_wt"] = pd.to_numeric(func[19], errors="coerce")
    func["bhatt_hetero_functional_class"] = func[20].astype("string")
    func["bhatt_wt_cotransfection_rescue"] = func[21].astype("string")

    cols_rel = [
        "protein_change",
        "bhatt_current_rel_wt",
        "bhatt_delta_vhalf",
        "bhatt_k_ratio",
        "bhatt_tau_ratio",
        "bhatt_phenotype",
    ]
    cols_func = [
        "protein_change",
        "bhatt_homo_current_rel_wt",
        "bhatt_homo_functional_class",
        "bhatt_hetero_current_rel_wt",
        "bhatt_hetero_functional_class",
        "bhatt_wt_cotransfection_rescue",
    ]

    het_retig_var = pd.read_excel(BENCHMARK_DIR / "bhatt_kcnq2_supp_table6.xlsx", sheet_name="Relative to variant control", header=None)
    het_retig_var = het_retig_var.iloc[16:].copy()
    het_retig_var = het_retig_var[het_retig_var[3].notna()].copy()
    het_retig_var["protein_change"] = het_retig_var[3].astype(str)
    het_retig_var["bhatt_retig_current_fold_vs_variant_ctrl"] = pd.to_numeric(het_retig_var[8], errors="coerce")
    het_retig_var["bhatt_retig_delta_vhalf_vs_variant_ctrl"] = pd.to_numeric(het_retig_var[14], errors="coerce")
    het_retig_var["bhatt_retig_k_ratio_vs_variant_ctrl"] = pd.to_numeric(het_retig_var[18], errors="coerce")
    het_retig_var["bhatt_retig_tau_ratio_vs_variant_ctrl"] = pd.to_numeric(het_retig_var[24], errors="coerce")

    het_retig_wt = pd.read_excel(BENCHMARK_DIR / "bhatt_kcnq2_supp_table6.xlsx", sheet_name="Relative to WT control", header=None)
    het_retig_wt = het_retig_wt.iloc[16:].copy()
    het_retig_wt = het_retig_wt[het_retig_wt[3].notna()].copy()
    het_retig_wt["protein_change"] = het_retig_wt[3].astype(str)
    het_retig_wt["bhatt_retig_current_rel_wt_ctrl"] = pd.to_numeric(het_retig_wt[8], errors="coerce")
    het_retig_wt["bhatt_retig_delta_vhalf_vs_wt_ctrl"] = pd.to_numeric(het_retig_wt[14], errors="coerce")
    het_retig_wt["bhatt_retig_k_ratio_vs_wt_ctrl"] = pd.to_numeric(het_retig_wt[18], errors="coerce")
    het_retig_wt["bhatt_retig_reaches_wt_current"] = het_retig_wt["bhatt_retig_current_rel_wt_ctrl"] >= 1.0

    return (
        rel[cols_rel]
        .merge(func[cols_func], on="protein_change", how="left")
        .merge(
            het_retig_var[
                [
                    "protein_change",
                    "bhatt_retig_current_fold_vs_variant_ctrl",
                    "bhatt_retig_delta_vhalf_vs_variant_ctrl",
                    "bhatt_retig_k_ratio_vs_variant_ctrl",
                    "bhatt_retig_tau_ratio_vs_variant_ctrl",
                ]
            ],
            on="protein_change",
            how="left",
        )
        .merge(
            het_retig_wt[
                [
                    "protein_change",
                    "bhatt_retig_current_rel_wt_ctrl",
                    "bhatt_retig_delta_vhalf_vs_wt_ctrl",
                    "bhatt_retig_k_ratio_vs_wt_ctrl",
                    "bhatt_retig_reaches_wt_current",
                ]
            ],
            on="protein_change",
            how="left",
        )
        .drop_duplicates("protein_change")
    )


def parse_vanoye_kcnq1() -> pd.DataFrame:
    path = BENCHMARK_DIR / "vanoye_kcnq1_dataset.xlsx"
    frames = []
    for sheet in ["Table S4a", "Table S4b"]:
        df = pd.read_excel(path, sheet_name=sheet, header=None)
        df = df.iloc[5:].copy()
        df = df[df[0].notna()].copy()
        df = df[~df[0].isin(["WT", "CHO"])].copy()
        parsed = pd.DataFrame(
            {
                "protein_change": df[0].astype(str),
                "vanoye_source_sheet": sheet,
                "vanoye_current_rel_wt": pd.to_numeric(df[5], errors="coerce"),
                "vanoye_delta_vhalf": pd.to_numeric(df[13], errors="coerce"),
                "vanoye_k": pd.to_numeric(df[21], errors="coerce"),
                "vanoye_tau_ms": pd.to_numeric(df[32], errors="coerce"),
            }
        )
        frames.append(parsed)
    return pd.concat(frames, ignore_index=True).drop_duplicates("protein_change")


def parse_brewer_kcnq1_mave() -> pd.DataFrame:
    path = BENCHMARK_DIR / "kcnq1_mave_scores.csv"
    if not path.exists():
        download(
            "https://raw.githubusercontent.com/GlazerLab/KCNQ1_MAVE/main/Analysis/KCNQ1Scores2_9.11.25.csv",
            path,
        )
    df = pd.read_csv(path)
    keep = [
        "mutation",
        "trafficking_score",
        "function_score",
        "het_trafficking_score",
        "het_function_score",
        "cluster_n6_name_FINAL",
        "clinvar_actual",
        "peakCurrent_lit",
        "deltaV12act_lit",
        "het_PeakCurrent_lit",
        "region",
    ]
    df = df[keep].rename(
        columns={
            "mutation": "protein_change",
            "trafficking_score": "brewer_trafficking_score",
            "function_score": "brewer_function_score",
            "het_trafficking_score": "brewer_het_trafficking_score",
            "het_function_score": "brewer_het_function_score",
            "cluster_n6_name_FINAL": "brewer_cluster",
            "clinvar_actual": "brewer_clinvar_actual",
            "peakCurrent_lit": "brewer_peak_current_lit",
            "deltaV12act_lit": "brewer_delta_vhalf_lit",
            "het_PeakCurrent_lit": "brewer_het_peak_current_lit",
            "region": "brewer_region",
        }
    )
    return df.drop_duplicates("protein_change").reset_index(drop=True)


def build_kcnq1_candidate_validation(df_top10: pd.DataFrame, vanoye_df: pd.DataFrame, brewer_df: pd.DataFrame) -> pd.DataFrame:
    candidates = df_top10[df_top10["gene"] == "KCNQ1"].copy()
    if candidates.empty:
        return pd.DataFrame()
    out = candidates.merge(vanoye_df, on="protein_change", how="left").merge(brewer_df, on="protein_change", how="left")

    def summarize(row: pd.Series) -> str:
        notes = []
        if pd.notna(row.get("vanoye_current_rel_wt")):
            notes.append("Vanoye overlap")
        else:
            notes.append("No Vanoye overlap")
        het_function = row.get("brewer_het_function_score")
        if pd.notna(het_function):
            if float(het_function) >= 0.75:
                notes.append("Brewer MAVE: retained heterozygous function")
            elif float(het_function) <= 0.25:
                notes.append("Brewer MAVE: severe heterozygous dysfunction")
            else:
                notes.append("Brewer MAVE: partial heterozygous dysfunction")
        else:
            notes.append("No Brewer MAVE match")
        return "; ".join(notes)

    out["validation_interpretation"] = out.apply(summarize, axis=1)
    return out[
        [
            "gene",
            "protein_change",
            "hgvs_c",
            "region_label",
            "structural_opportunity_score",
            "tractability_modifier",
            "rescue_priority_score",
            "vanoye_source_sheet",
            "vanoye_current_rel_wt",
            "vanoye_delta_vhalf",
            "vanoye_k",
            "vanoye_tau_ms",
            "brewer_trafficking_score",
            "brewer_function_score",
            "brewer_het_trafficking_score",
            "brewer_het_function_score",
            "brewer_cluster",
            "brewer_clinvar_actual",
            "validation_interpretation",
        ]
    ].sort_values("rescue_priority_score", ascending=False)


def build_ml277_status_table() -> pd.DataFrame:
    return pd.DataFrame(
        [
            {
                "agent": "ML277",
                "reported_status": "No registered ML277 study page identified via ClinicalTrials.gov search term=ML277 as of 2026-03-22",
                "evidence_type": "ClinicalTrials.gov search check",
                "search_url": "https://clinicaltrials.gov/search?term=ML277",
                "landscape_note": "Current KCNQ1 rescue literature includes preclinical trafficking corrector VU0494372 rather than an ML277 clinical trial",
                "supporting_source": "JCI Insight 2026 article 201297; ClinicalTrials.gov search",
            }
        ]
    )


def benchmark_overlap(df_clinical: pd.DataFrame, benchmark_df: pd.DataFrame, gene: str, benchmark_prefix: str) -> tuple[pd.DataFrame, pd.DataFrame]:
    gene_df = df_clinical[df_clinical["gene"] == gene].copy()
    existing_benchmark_cols = [col for col in gene_df.columns if col.startswith(f"{benchmark_prefix}_")]
    if existing_benchmark_cols:
        gene_df = gene_df.drop(columns=existing_benchmark_cols)
    overlap = gene_df.merge(benchmark_df, on="protein_change", how="inner")
    rows = [
        {
            "gene": gene,
            "benchmark": benchmark_prefix,
            "clinical_variants": int(len(gene_df)),
            "benchmark_variants": int(len(benchmark_df)),
            "overlap_variants": int(len(overlap)),
        }
    ]
    if len(overlap) >= 3:
        value_col = f"{benchmark_prefix}_current_rel_wt"
        rho, p_value = stats.spearmanr(overlap["rescue_score"], overlap[value_col], nan_policy="omit")
        rows[0]["spearman_rho_rescue_vs_current_rel_wt"] = rho
        rows[0]["spearman_p_value"] = p_value
    else:
        rows[0]["spearman_rho_rescue_vs_current_rel_wt"] = np.nan
        rows[0]["spearman_p_value"] = np.nan
    return overlap, pd.DataFrame(rows)


def annotate_bhatt_mechanism(overlap: pd.DataFrame) -> pd.DataFrame:
    out = overlap.copy()
    mechanism_map = {
        "-D N-": "dominant_negative",
        "-f R-": "full_rescue",
        "-p R-": "partial_rescue",
        "none": "no_rescue_annotated",
    }
    out["bhatt_mechanism_class"] = out["bhatt_wt_cotransfection_rescue"].map(mechanism_map).fillna("unassigned")
    out["bhatt_retig_responsive_current50"] = out["bhatt_retig_current_fold_vs_variant_ctrl"] >= 1.5
    out["bhatt_retig_responsive_vhalf10"] = out["bhatt_retig_delta_vhalf_vs_variant_ctrl"].abs() >= 10.0
    out["bhatt_retig_responsive_any"] = out["bhatt_retig_responsive_current50"].fillna(False) | out[
        "bhatt_retig_responsive_vhalf10"
    ].fillna(False)
    return out


def _score_endpoint_auc(overlap: pd.DataFrame, score_col: str, endpoint_col: str) -> float:
    df = overlap.dropna(subset=[score_col, endpoint_col]).copy()
    if len(df) < 3 or len(df[endpoint_col].astype(int).unique()) < 2:
        return np.nan
    return float(roc_auc_score(df[endpoint_col].astype(int), df[score_col]))


def _score_endpoint_pr_auc(overlap: pd.DataFrame, score_col: str, endpoint_col: str) -> float:
    df = overlap.dropna(subset=[score_col, endpoint_col]).copy()
    if len(df) < 3 or len(df[endpoint_col].astype(int).unique()) < 2:
        return np.nan
    return float(average_precision_score(df[endpoint_col].astype(int), df[score_col]))


def _bootstrap_auc_ci(
    overlap: pd.DataFrame,
    score_col: str,
    endpoint_col: str,
    n_boot: int = 10_000,
    seed: int = 42,
) -> tuple[float, float]:
    df = overlap.dropna(subset=[score_col, endpoint_col]).copy()
    if len(df) < 8 or len(df[endpoint_col].astype(int).unique()) < 2:
        return np.nan, np.nan
    rng = np.random.default_rng(seed)
    y = df[endpoint_col].astype(int).to_numpy()
    scores = df[score_col].to_numpy(dtype=float)
    aucs = []
    idx = np.arange(len(df))
    for _ in range(n_boot):
        sample_idx = rng.choice(idx, size=len(idx), replace=True)
        y_boot = y[sample_idx]
        if len(np.unique(y_boot)) < 2:
            continue
        aucs.append(roc_auc_score(y_boot, scores[sample_idx]))
    if len(aucs) < 100:
        return np.nan, np.nan
    lo, hi = np.percentile(aucs, [2.5, 97.5])
    return float(lo), float(hi)


def _score_endpoint_spearman(overlap: pd.DataFrame, score_col: str, endpoint_col: str) -> tuple[float, float]:
    df = overlap.dropna(subset=[score_col, endpoint_col]).copy()
    if len(df) < 3:
        return np.nan, np.nan
    rho, p_value = stats.spearmanr(df[score_col], df[endpoint_col], nan_policy="omit")
    return float(rho), float(p_value)


def _compute_midrank(x: np.ndarray) -> np.ndarray:
    order = np.argsort(x)
    sorted_x = x[order]
    n = len(x)
    midranks = np.zeros(n, dtype=float)
    i = 0
    while i < n:
        j = i
        while j < n and sorted_x[j] == sorted_x[i]:
            j += 1
        midranks[i:j] = 0.5 * (i + j - 1) + 1.0
        i = j
    out = np.empty(n, dtype=float)
    out[order] = midranks
    return out


def _fast_delong(predictions_sorted_transposed: np.ndarray, label_1_count: int) -> tuple[np.ndarray, np.ndarray]:
    m = label_1_count
    n = predictions_sorted_transposed.shape[1] - m
    positive_examples = predictions_sorted_transposed[:, :m]
    negative_examples = predictions_sorted_transposed[:, m:]
    k = predictions_sorted_transposed.shape[0]

    tx = np.empty((k, m), dtype=float)
    ty = np.empty((k, n), dtype=float)
    tz = np.empty((k, m + n), dtype=float)
    for r in range(k):
        tx[r, :] = _compute_midrank(positive_examples[r, :])
        ty[r, :] = _compute_midrank(negative_examples[r, :])
        tz[r, :] = _compute_midrank(predictions_sorted_transposed[r, :])

    aucs = tz[:, :m].sum(axis=1) / (m * n) - (m + 1.0) / (2.0 * n)
    v01 = (tz[:, :m] - tx) / n
    v10 = 1.0 - (tz[:, m:] - ty) / m
    sx = np.cov(v01)
    sy = np.cov(v10)
    delong_cov = sx / m + sy / n
    return aucs, delong_cov


def _delong_roc_test(y_true: np.ndarray, pred_a: np.ndarray, pred_b: np.ndarray) -> tuple[float, float, float, float]:
    y_true = np.asarray(y_true, dtype=int)
    pred_a = np.asarray(pred_a, dtype=float)
    pred_b = np.asarray(pred_b, dtype=float)
    if len(np.unique(y_true)) < 2:
        return np.nan, np.nan, np.nan, np.nan
    order = np.argsort(-y_true)
    label_1_count = int(y_true.sum())
    preds = np.vstack([pred_a, pred_b])[:, order]
    aucs, cov = _fast_delong(preds, label_1_count)
    diff = float(aucs[0] - aucs[1])
    var = float(cov[0, 0] + cov[1, 1] - 2.0 * cov[0, 1])
    if var <= 0:
        return float(aucs[0]), float(aucs[1]), diff, np.nan
    z = abs(diff) / np.sqrt(var)
    p_value = 2.0 * stats.norm.sf(z)
    return float(aucs[0]), float(aucs[1]), diff, float(p_value)


def _classification_metrics(y_true: np.ndarray, y_pred: np.ndarray) -> dict[str, float]:
    y_true = np.asarray(y_true, dtype=int)
    y_pred = np.asarray(y_pred, dtype=int)
    tp = int(((y_true == 1) & (y_pred == 1)).sum())
    tn = int(((y_true == 0) & (y_pred == 0)).sum())
    fp = int(((y_true == 0) & (y_pred == 1)).sum())
    fn = int(((y_true == 1) & (y_pred == 0)).sum())

    sensitivity = tp / (tp + fn) if (tp + fn) else np.nan
    specificity = tn / (tn + fp) if (tn + fp) else np.nan
    ppv = tp / (tp + fp) if (tp + fp) else np.nan
    npv = tn / (tn + fn) if (tn + fn) else np.nan
    accuracy = (tp + tn) / len(y_true) if len(y_true) else np.nan
    f1 = 2 * tp / (2 * tp + fp + fn) if (2 * tp + fp + fn) else np.nan
    balanced_accuracy = np.nanmean([sensitivity, specificity])
    denom = np.sqrt((tp + fp) * (tp + fn) * (tn + fp) * (tn + fn))
    mcc = ((tp * tn) - (fp * fn)) / denom if denom else np.nan
    prevalence = float(y_true.mean()) if len(y_true) else np.nan
    return {
        "tp": tp,
        "tn": tn,
        "fp": fp,
        "fn": fn,
        "sensitivity": sensitivity,
        "specificity": specificity,
        "ppv": ppv,
        "npv": npv,
        "accuracy": accuracy,
        "f1_score": f1,
        "balanced_accuracy": balanced_accuracy,
        "mcc": mcc,
        "prevalence": prevalence,
    }


def benchmark_bhatt_retigabine(df_clinical: pd.DataFrame, bhatt_df: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    gene_df = df_clinical[df_clinical["gene"] == "KCNQ2"].copy()
    existing_bhatt_cols = [col for col in gene_df.columns if col.startswith("bhatt_")]
    if existing_bhatt_cols:
        gene_df = gene_df.drop(columns=existing_bhatt_cols)
    overlap = gene_df.merge(bhatt_df, on="protein_change", how="inner")
    overlap = annotate_bhatt_mechanism(overlap)
    summary = {
        "gene": "KCNQ2",
        "benchmark": "bhatt_retigabine",
        "overlap_variants": int(len(overlap)),
        "n_with_retig_current_rel_wt": int(overlap["bhatt_retig_current_rel_wt_ctrl"].notna().sum()),
        "n_reaches_wt_current": int(overlap["bhatt_retig_reaches_wt_current"].fillna(False).sum()),
        "n_responsive_current50": int(overlap["bhatt_retig_responsive_current50"].fillna(False).sum()),
        "n_responsive_vhalf10": int(overlap["bhatt_retig_responsive_vhalf10"].fillna(False).sum()),
        "n_responsive_any": int(overlap["bhatt_retig_responsive_any"].fillna(False).sum()),
    }
    score_cols = ["rescue_score"]
    if "rescue_priority_score" in overlap.columns:
        score_cols.append("rescue_priority_score")
    for score_col in score_cols:
        prefix = "priority" if score_col == "rescue_priority_score" else "structural"
        rho_current, p_current = _score_endpoint_spearman(overlap, score_col, "bhatt_retig_current_rel_wt_ctrl")
        summary[f"spearman_rho_{prefix}_vs_retig_current_rel_wt"] = rho_current
        summary[f"spearman_p_{prefix}_vs_retig_current_rel_wt"] = p_current
        rho_vhalf, p_vhalf = _score_endpoint_spearman(overlap, score_col, "bhatt_retig_delta_vhalf_vs_variant_ctrl")
        summary[f"spearman_rho_{prefix}_vs_retig_delta_vhalf_variant_ctrl"] = rho_vhalf
        summary[f"spearman_p_{prefix}_vs_retig_delta_vhalf_variant_ctrl"] = p_vhalf
        summary[f"auc_{prefix}_predict_reaches_wt_current"] = _score_endpoint_auc(
            overlap, score_col, "bhatt_retig_reaches_wt_current"
        )
        summary[f"auc_{prefix}_predict_responsive_current50"] = _score_endpoint_auc(
            overlap, score_col, "bhatt_retig_responsive_current50"
        )
        summary[f"auc_{prefix}_predict_responsive_any"] = _score_endpoint_auc(
            overlap, score_col, "bhatt_retig_responsive_any"
        )
        summary[f"pr_auc_{prefix}_predict_reaches_wt_current"] = _score_endpoint_pr_auc(
            overlap, score_col, "bhatt_retig_reaches_wt_current"
        )
        summary[f"pr_auc_{prefix}_predict_responsive_current50"] = _score_endpoint_pr_auc(
            overlap, score_col, "bhatt_retig_responsive_current50"
        )
        summary[f"pr_auc_{prefix}_predict_responsive_any"] = _score_endpoint_pr_auc(
            overlap, score_col, "bhatt_retig_responsive_any"
        )
        ci_lo, ci_hi = _bootstrap_auc_ci(overlap, score_col, "bhatt_retig_reaches_wt_current")
        summary[f"auc_{prefix}_predict_reaches_wt_current_ci_lo"] = ci_lo
        summary[f"auc_{prefix}_predict_reaches_wt_current_ci_hi"] = ci_hi
        ci_lo, ci_hi = _bootstrap_auc_ci(overlap, score_col, "bhatt_retig_responsive_current50")
        summary[f"auc_{prefix}_predict_responsive_current50_ci_lo"] = ci_lo
        summary[f"auc_{prefix}_predict_responsive_current50_ci_hi"] = ci_hi
        ci_lo, ci_hi = _bootstrap_auc_ci(overlap, score_col, "bhatt_retig_responsive_any")
        summary[f"auc_{prefix}_predict_responsive_any_ci_lo"] = ci_lo
        summary[f"auc_{prefix}_predict_responsive_any_ci_hi"] = ci_hi
    return overlap, pd.DataFrame([summary])


def benchmark_bhatt_operating_point(
    df_clinical: pd.DataFrame,
    bhatt_df: pd.DataFrame,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    gene_df = df_clinical[df_clinical["gene"] == "KCNQ2"].copy()
    existing_bhatt_cols = [col for col in gene_df.columns if col.startswith("bhatt_")]
    if existing_bhatt_cols:
        gene_df = gene_df.drop(columns=existing_bhatt_cols)
    overlap = gene_df.merge(bhatt_df, on="protein_change", how="inner")
    overlap = annotate_bhatt_mechanism(overlap)

    endpoint_col = "bhatt_retig_responsive_any"
    endpoint_label = "responsive_any"
    model_specs = [
        {
            "model_name": "rescue_priority_score",
            "score_col": "rescue_priority_score",
            "threshold": 0.50,
            "threshold_rule": "pre-specified fixed priority threshold",
        },
        {
            "model_name": "alphamissense_score",
            "score_col": "alphamissense_score",
            "threshold": 0.75,
            "threshold_rule": "standard AlphaMissense pathogenicity threshold",
        },
    ]

    metrics_rows = []
    confusion_rows = []
    for spec in model_specs:
        score_col = spec["score_col"]
        if score_col not in overlap.columns:
            continue
        use = overlap.dropna(subset=[score_col, endpoint_col]).copy()
        if use.empty or len(use[endpoint_col].astype(int).unique()) < 2:
            continue
        y_true = use[endpoint_col].astype(int).to_numpy()
        y_pred = (use[score_col].to_numpy(dtype=float) >= spec["threshold"]).astype(int)
        metrics = _classification_metrics(y_true, y_pred)
        metrics_rows.append(
            {
                "endpoint": endpoint_label,
                "endpoint_column": endpoint_col,
                "model_name": spec["model_name"],
                "score_column": score_col,
                "threshold": spec["threshold"],
                "threshold_rule": spec["threshold_rule"],
                "n_variants": int(len(use)),
                "n_positive": int(y_true.sum()),
                "n_negative": int((1 - y_true).sum()),
                **metrics,
            }
        )
        confusion_rows.extend(
            [
                {
                    "endpoint": endpoint_label,
                    "model_name": spec["model_name"],
                    "threshold": spec["threshold"],
                    "actual_class": "positive",
                    "predicted_class": "positive",
                    "count": metrics["tp"],
                },
                {
                    "endpoint": endpoint_label,
                    "model_name": spec["model_name"],
                    "threshold": spec["threshold"],
                    "actual_class": "positive",
                    "predicted_class": "negative",
                    "count": metrics["fn"],
                },
                {
                    "endpoint": endpoint_label,
                    "model_name": spec["model_name"],
                    "threshold": spec["threshold"],
                    "actual_class": "negative",
                    "predicted_class": "positive",
                    "count": metrics["fp"],
                },
                {
                    "endpoint": endpoint_label,
                    "model_name": spec["model_name"],
                    "threshold": spec["threshold"],
                    "actual_class": "negative",
                    "predicted_class": "negative",
                    "count": metrics["tn"],
                },
            ]
        )
    return pd.DataFrame(metrics_rows), pd.DataFrame(confusion_rows)


def benchmark_bhatt_mechanism_strata(overlap: pd.DataFrame) -> pd.DataFrame:
    rows = []
    score_cols = ["rescue_score"]
    if "rescue_priority_score" in overlap.columns:
        score_cols.append("rescue_priority_score")
    for mechanism, group in overlap.groupby("bhatt_mechanism_class", dropna=False):
        row = {
            "mechanism_class": mechanism,
            "n_variants": int(len(group)),
            "median_structural_opportunity_score": float(group["rescue_score"].median()),
            "median_tractability_modifier": float(group["tractability_modifier"].median()) if "tractability_modifier" in group.columns else np.nan,
            "median_rescue_priority_score": float(group["rescue_priority_score"].median()) if "rescue_priority_score" in group.columns else np.nan,
            "median_retig_current_rel_wt_ctrl": float(group["bhatt_retig_current_rel_wt_ctrl"].median()),
            "median_retig_current_fold_vs_variant_ctrl": float(group["bhatt_retig_current_fold_vs_variant_ctrl"].median()),
            "responsive_current50_rate": float(group["bhatt_retig_responsive_current50"].fillna(False).mean()),
            "responsive_any_rate": float(group["bhatt_retig_responsive_any"].fillna(False).mean()),
            "reaches_wt_current_rate": float(group["bhatt_retig_reaches_wt_current"].fillna(False).mean()),
        }
        for score_col in score_cols:
            prefix = "priority" if score_col == "rescue_priority_score" else "structural"
            rho, p_value = _score_endpoint_spearman(group, score_col, "bhatt_retig_current_rel_wt_ctrl")
            row[f"spearman_rho_{prefix}_vs_retig_current_rel_wt"] = rho
            row[f"spearman_p_{prefix}_vs_retig_current_rel_wt"] = p_value
        rows.append(row)
    return pd.DataFrame(rows).sort_values(["n_variants", "mechanism_class"], ascending=[False, True]).reset_index(drop=True)


def compute_alphamissense_concordance(df_clinical: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    use = df_clinical.dropna(subset=["alphamissense_score", "structural_opportunity_score"]).copy()
    if use.empty:
        return pd.DataFrame(), pd.DataFrame()
    use["variant_label"] = use["gene"] + " " + use["protein_change"]
    use["clinvar_group"] = np.where(
        use["clinical_significance"].isin(["Pathogenic", "Likely pathogenic", "Pathogenic/Likely pathogenic"]),
        "P/LP",
        np.where(
            use["clinical_significance"].isin(["Uncertain significance", "Conflicting classifications of pathogenicity"]),
            "VUS/conflicting",
            "Other",
        ),
    )
    am_cutoff = 0.75
    structural_cutoff = 0.70
    use["quadrant"] = np.select(
        [
            (use["alphamissense_score"] >= am_cutoff) & (use["structural_opportunity_score"] >= structural_cutoff),
            (use["alphamissense_score"] < am_cutoff) & (use["structural_opportunity_score"] >= structural_cutoff),
            (use["alphamissense_score"] >= am_cutoff) & (use["structural_opportunity_score"] < structural_cutoff),
        ],
        [
            "pathogenic_and_druggable",
            "structurally_druggable_but_lower_am",
            "pathogenic_but_structurally_peripheral",
        ],
        default="background_low_low",
    )
    spearman_rho, spearman_p = stats.spearmanr(use["alphamissense_score"], use["structural_opportunity_score"], nan_policy="omit")
    pearson_r, pearson_p = stats.pearsonr(use["alphamissense_score"], use["structural_opportunity_score"])
    summary = pd.DataFrame(
        [
            {
                "n_variants": int(len(use)),
                "alphamissense_high_cutoff": am_cutoff,
                "structural_high_cutoff": structural_cutoff,
                "spearman_rho": float(spearman_rho),
                "spearman_p_value": float(spearman_p),
                "pearson_r": float(pearson_r),
                "pearson_p_value": float(pearson_p),
                "n_pathogenic_and_druggable": int((use["quadrant"] == "pathogenic_and_druggable").sum()),
                "n_structurally_druggable_but_lower_am": int((use["quadrant"] == "structurally_druggable_but_lower_am").sum()),
                "n_pathogenic_but_structurally_peripheral": int((use["quadrant"] == "pathogenic_but_structurally_peripheral").sum()),
                "n_background_low_low": int((use["quadrant"] == "background_low_low").sum()),
            }
        ]
    )
    return use, summary


def benchmark_alphamissense_vs_structural(df_clinical: pd.DataFrame, bhatt_df: pd.DataFrame) -> pd.DataFrame:
    overlap = df_clinical[df_clinical["gene"] == "KCNQ2"].merge(
        annotate_bhatt_mechanism(bhatt_df.copy()),
        on="protein_change",
        how="inner",
    )
    endpoints = {
        "responsive_current50": "bhatt_retig_responsive_current50",
        "responsive_any": "bhatt_retig_responsive_any",
        "reaches_wt_current": "bhatt_retig_reaches_wt_current",
    }
    rows = []
    for endpoint_name, endpoint_col in endpoints.items():
        use = overlap.dropna(subset=["alphamissense_score", "rescue_priority_score", endpoint_col]).copy()
        if len(use) < 3 or len(use[endpoint_col].astype(int).unique()) < 2:
            continue
        z = StandardScaler().fit_transform(use[["alphamissense_score", "rescue_priority_score"]])
        use["combined_mean_z"] = z.mean(axis=1)
        models = {
            "alphamissense_only": "alphamissense_score",
            "rescue_priority_only": "rescue_priority_score",
            "combined_mean_z": "combined_mean_z",
        }
        for model_name, score_col in models.items():
            rows.append(
                {
                    "model_name": model_name,
                    "score_column": score_col,
                    "endpoint": endpoint_name,
                    "n_variants_used": int(len(use)),
                    "auc": float(roc_auc_score(use[endpoint_col].astype(int), use[score_col])),
                }
            )
    return pd.DataFrame(rows).sort_values(["endpoint", "auc"], ascending=[True, False]).reset_index(drop=True)


def benchmark_bhatt_delong(df_clinical: pd.DataFrame, bhatt_df: pd.DataFrame) -> pd.DataFrame:
    overlap = df_clinical[df_clinical["gene"] == "KCNQ2"].merge(
        annotate_bhatt_mechanism(bhatt_df.copy()),
        on="protein_change",
        how="inner",
    )
    endpoints = {
        "responsive_any": "bhatt_retig_responsive_any",
        "responsive_current50": "bhatt_retig_responsive_current50",
        "reaches_wt_current": "bhatt_retig_reaches_wt_current",
    }
    comparison_specs = [
        ("rescue_priority_vs_alphamissense", "rescue_priority_score", "alphamissense_score"),
        ("structural_vs_alphamissense", "structural_opportunity_score", "alphamissense_score"),
    ]
    rows = []
    for endpoint_name, endpoint_col in endpoints.items():
        for comparison_name, model_a, model_b in comparison_specs:
            use = overlap.dropna(subset=[model_a, model_b, endpoint_col]).copy()
            if len(use) < 8 or len(use[endpoint_col].astype(int).unique()) < 2:
                continue
            auc_a, auc_b, auc_diff, p_value = _delong_roc_test(
                use[endpoint_col].astype(int).to_numpy(),
                use[model_a].to_numpy(dtype=float),
                use[model_b].to_numpy(dtype=float),
            )
            rows.append(
                {
                    "endpoint": endpoint_name,
                    "endpoint_column": endpoint_col,
                    "comparison": comparison_name,
                    "model_a": model_a,
                    "model_b": model_b,
                    "n_variants_used": int(len(use)),
                    "auc_model_a": auc_a,
                    "auc_model_b": auc_b,
                    "auc_difference_a_minus_b": auc_diff,
                    "delong_p_value": p_value,
                    "better_model": model_a if auc_diff > 0 else model_b,
                }
            )
    return pd.DataFrame(rows).sort_values(["endpoint", "comparison"]).reset_index(drop=True)


def benchmark_alphamissense_high_subsets(df_clinical: pd.DataFrame, bhatt_df: pd.DataFrame) -> pd.DataFrame:
    overlap = df_clinical[df_clinical["gene"] == "KCNQ2"].merge(
        annotate_bhatt_mechanism(bhatt_df.copy()),
        on="protein_change",
        how="inner",
    )
    rows = []
    for cutoff in [0.90, 0.95, 0.98, 0.99]:
        subset = overlap[overlap["alphamissense_score"] >= cutoff].copy()
        for endpoint_name, endpoint_col in {
            "responsive_any": "bhatt_retig_responsive_any",
            "responsive_current50": "bhatt_retig_responsive_current50",
            "reaches_wt_current": "bhatt_retig_reaches_wt_current",
        }.items():
            use = subset.dropna(subset=["alphamissense_score", "structural_opportunity_score", "rescue_priority_score", endpoint_col]).copy()
            if len(use) < 8 or len(use[endpoint_col].astype(int).unique()) < 2:
                continue
            for model_name, score_col in {
                "alphamissense_only": "alphamissense_score",
                "structural_only": "structural_opportunity_score",
                "rescue_priority_only": "rescue_priority_score",
            }.items():
                rows.append(
                    {
                        "alphamissense_cutoff": cutoff,
                        "endpoint": endpoint_name,
                        "n_variants_used": int(len(use)),
                        "model_name": model_name,
                        "auc": float(roc_auc_score(use[endpoint_col].astype(int), use[score_col])),
                    }
                )
    return pd.DataFrame(rows).sort_values(["alphamissense_cutoff", "endpoint", "auc"], ascending=[True, True, False]).reset_index(drop=True)


def benchmark_revel_vs_structural(df_clinical: pd.DataFrame, bhatt_df: pd.DataFrame, revel_df: pd.DataFrame) -> pd.DataFrame:
    overlap = df_clinical[df_clinical["gene"] == "KCNQ2"].merge(
        annotate_bhatt_mechanism(bhatt_df.copy()),
        on="protein_change",
        how="inner",
    )
    overlap = overlap.merge(revel_df[["protein_change", "revel_score"]], on="protein_change", how="left")
    endpoints = {
        "responsive_current50": "bhatt_retig_responsive_current50",
        "responsive_any": "bhatt_retig_responsive_any",
        "reaches_wt_current": "bhatt_retig_reaches_wt_current",
    }
    rows = []
    for endpoint_name, endpoint_col in endpoints.items():
        use = overlap.dropna(subset=["revel_score", "rescue_priority_score", endpoint_col]).copy()
        if len(use) < 3 or len(use[endpoint_col].astype(int).unique()) < 2:
            continue
        models = {
            "revel_only": "revel_score",
            "alphamissense_only": "alphamissense_score",
            "rescue_priority_only": "rescue_priority_score",
            "structural_only": "structural_opportunity_score",
        }
        for model_name, score_col in models.items():
            if score_col not in use.columns:
                continue
            sub = use.dropna(subset=[score_col]).copy()
            if len(sub) < 3 or len(sub[endpoint_col].astype(int).unique()) < 2:
                continue
            rows.append(
                {
                    "endpoint": endpoint_name,
                    "model_name": model_name,
                    "score_column": score_col,
                    "n_variants_used": int(len(sub)),
                    "auc": float(roc_auc_score(sub[endpoint_col].astype(int), sub[score_col])),
                }
            )
    return pd.DataFrame(rows).sort_values(["endpoint", "auc"], ascending=[True, False]).reset_index(drop=True)


def _assign_brewer_mechanism_like(row: pd.Series) -> str:
    trafficking = row.get("brewer_trafficking_score", row.get("trafficking_score"))
    het_trafficking = row.get("brewer_het_trafficking_score", row.get("het_trafficking_score"))
    peak_current = row.get("brewer_peak_current_lit", row.get("peakCurrent_lit"))
    het_peak_current = row.get("brewer_het_peak_current_lit", row.get("het_PeakCurrent_lit"))
    delta_vhalf = row.get("brewer_delta_vhalf_lit", row.get("deltaV12act_lit"))
    function_score = row.get("brewer_function_score", row.get("function_score"))
    het_function_score = row.get("brewer_het_function_score", row.get("het_function_score"))
    if pd.notna(delta_vhalf) and abs(float(delta_vhalf)) >= 10 and (pd.isna(trafficking) or float(trafficking) >= 0.5):
        return "gating_shift"
    low_trafficking = (pd.notna(trafficking) and float(trafficking) < 0.5) or (
        pd.notna(het_trafficking) and float(het_trafficking) < 0.5
    )
    low_peak = (pd.notna(peak_current) and float(peak_current) < 0.5) or (
        pd.notna(het_peak_current) and float(het_peak_current) < 0.5
    )
    low_function = (pd.notna(function_score) and float(function_score) < 0.5) or (
        pd.notna(het_function_score) and float(het_function_score) < 0.5
    )
    if low_trafficking and not (pd.notna(delta_vhalf) and abs(float(delta_vhalf)) >= 10):
        return "instability_mistrafficking"
    if low_peak or low_function:
        return "conductance_loss"
    return "other_or_mixed"


def benchmark_brewer_kcnq1_mechanisms(df_clinical: pd.DataFrame, brewer_df: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    kcnq1 = df_clinical[df_clinical["gene"] == "KCNQ1"].copy()
    overlap = kcnq1.merge(brewer_df, on="protein_change", how="inner")
    overlap["brewer_mechanism_like"] = overlap.apply(_assign_brewer_mechanism_like, axis=1)
    lit_overlap = overlap[
        overlap[["brewer_peak_current_lit", "brewer_delta_vhalf_lit", "brewer_het_peak_current_lit"]].notna().any(axis=1)
    ].copy()

    summary_rows = []
    for mechanism, group in lit_overlap.groupby("brewer_mechanism_like", dropna=False):
        row = {
            "mechanism_like": mechanism,
            "n_variants": int(len(group)),
            "median_structural_opportunity_score": float(group["structural_opportunity_score"].median()),
            "median_rescue_priority_score": float(group["rescue_priority_score"].median()),
            "median_peak_current_lit": float(group["brewer_peak_current_lit"].median()) if group["brewer_peak_current_lit"].notna().any() else np.nan,
            "median_delta_vhalf_lit": float(group["brewer_delta_vhalf_lit"].median()) if group["brewer_delta_vhalf_lit"].notna().any() else np.nan,
            "median_het_peak_current_lit": float(group["brewer_het_peak_current_lit"].median()) if group["brewer_het_peak_current_lit"].notna().any() else np.nan,
            "median_trafficking_score": float(group["brewer_trafficking_score"].median()) if group["brewer_trafficking_score"].notna().any() else np.nan,
            "median_het_function_score": float(group["brewer_het_function_score"].median()) if group["brewer_het_function_score"].notna().any() else np.nan,
        }
        for score_col in ["structural_opportunity_score", "rescue_priority_score"]:
            sub = group.dropna(subset=[score_col, "brewer_het_function_score"])
            rho, p_value = _score_endpoint_spearman(sub, score_col, "brewer_het_function_score")
            prefix = "priority" if score_col == "rescue_priority_score" else "structural"
            row[f"spearman_rho_{prefix}_vs_het_function"] = rho
            row[f"spearman_p_{prefix}_vs_het_function"] = p_value
        summary_rows.append(row)
    summary = pd.DataFrame(summary_rows).sort_values("median_rescue_priority_score", ascending=False).reset_index(drop=True)

    omnibus_rows = []
    for score_col in ["structural_opportunity_score", "rescue_priority_score"]:
        groups = [g[score_col].dropna().to_numpy() for _, g in lit_overlap.groupby("brewer_mechanism_like") if len(g[score_col].dropna()) >= 3]
        if len(groups) >= 2:
            test = stats.kruskal(*groups)
            omnibus_rows.append({"score_col": score_col, "kruskal_statistic": float(test.statistic), "p_value": float(test.pvalue)})
    omnibus = pd.DataFrame(omnibus_rows)
    return overlap.reset_index(drop=True), summary, omnibus


def benchmark_bhatt_heteromeric_rescue(df_clinical: pd.DataFrame, bhatt_df: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    gene_df = df_clinical[df_clinical["gene"] == "KCNQ2"].copy()
    existing_bhatt_cols = [col for col in gene_df.columns if col.startswith("bhatt_")]
    if existing_bhatt_cols:
        gene_df = gene_df.drop(columns=existing_bhatt_cols)
    overlap = gene_df.merge(
        annotate_bhatt_mechanism(bhatt_df.copy()),
        on="protein_change",
        how="inner",
    )
    overlap["heteromeric_rescue_gain_over_homo"] = overlap["bhatt_retig_current_rel_wt_ctrl"] / overlap["bhatt_homo_current_rel_wt"]
    overlap["heteromeric_rescue_gain_over_hetero_baseline"] = overlap["bhatt_retig_current_rel_wt_ctrl"] / overlap["bhatt_hetero_current_rel_wt"]
    overlap = overlap.replace([np.inf, -np.inf], np.nan)

    summary_rows = []
    endpoint_specs = [
        ("bhatt_retig_current_rel_wt_ctrl", "hetero_retig_current_rel_wt"),
        ("bhatt_hetero_current_rel_wt", "hetero_baseline_current_rel_wt"),
        ("bhatt_homo_current_rel_wt", "homo_baseline_current_rel_wt"),
        ("bhatt_retig_current_fold_vs_variant_ctrl", "retig_fold_change_vs_variant"),
        ("heteromeric_rescue_gain_over_homo", "heteromeric_gain_over_homo"),
        ("heteromeric_rescue_gain_over_hetero_baseline", "heteromeric_gain_over_hetero_baseline"),
    ]
    for score_col in ["structural_opportunity_score", "rescue_priority_score"]:
        prefix = "priority" if score_col == "rescue_priority_score" else "structural"
        for endpoint_col, endpoint_label in endpoint_specs:
            rho, p_value = _score_endpoint_spearman(overlap, score_col, endpoint_col)
            sub = overlap.dropna(subset=[score_col, endpoint_col])
            summary_rows.append(
                {
                    "score_model": prefix,
                    "endpoint": endpoint_label,
                    "n_variants": int(len(sub)),
                    "spearman_rho": rho,
                    "spearman_p_value": p_value,
                    "median_endpoint_value": float(sub[endpoint_col].median()) if len(sub) else np.nan,
                }
            )

    mechanism_rows = []
    for mechanism, group in overlap.groupby("bhatt_mechanism_class", dropna=False):
        mechanism_rows.append(
            {
                "mechanism_class": mechanism,
                "n_variants": int(len(group)),
                "median_hetero_retig_current_rel_wt": float(group["bhatt_retig_current_rel_wt_ctrl"].median()) if group["bhatt_retig_current_rel_wt_ctrl"].notna().any() else np.nan,
                "median_heteromeric_gain_over_homo": float(group["heteromeric_rescue_gain_over_homo"].median()) if group["heteromeric_rescue_gain_over_homo"].notna().any() else np.nan,
                "median_heteromeric_gain_over_hetero_baseline": float(group["heteromeric_rescue_gain_over_hetero_baseline"].median()) if group["heteromeric_rescue_gain_over_hetero_baseline"].notna().any() else np.nan,
                "spearman_rho_priority_vs_hetero_retig_current_rel_wt": _score_endpoint_spearman(group, "rescue_priority_score", "bhatt_retig_current_rel_wt_ctrl")[0],
                "spearman_p_priority_vs_hetero_retig_current_rel_wt": _score_endpoint_spearman(group, "rescue_priority_score", "bhatt_retig_current_rel_wt_ctrl")[1],
            }
        )
    mechanism_summary = pd.DataFrame(mechanism_rows).sort_values("n_variants", ascending=False).reset_index(drop=True)
    return overlap.reset_index(drop=True), pd.DataFrame(summary_rows), mechanism_summary


def _get_uniprot_search_rows(query: str) -> pd.DataFrame:
    url = "https://rest.uniprot.org/uniprotkb/search"
    resp = session().get(
        url,
        params={
            "query": query,
            "fields": "accession,gene_primary,protein_name,organism_name,reviewed,length",
            "format": "tsv",
            "size": 10,
        },
        timeout=60,
    )
    resp.raise_for_status()
    lines = [line for line in resp.text.splitlines() if line.strip()]
    if len(lines) <= 1:
        return pd.DataFrame(columns=["Entry", "Gene Names (primary)", "Protein names", "Organism", "Reviewed", "Length"])
    data = [line.split("\t") for line in lines[1:]]
    return pd.DataFrame(data, columns=lines[0].split("\t"))


def _find_ortholog_accession(gene: str, species_label: str) -> tuple[str | None, str]:
    if species_label == "human":
        return UNIPROT_IDS[gene], "reviewed_human_reference"
    taxid = ORTHOLOG_SPECIES[species_label]
    queries = [
        f"gene_exact:{gene} AND organism_id:{taxid} AND reviewed:true",
        f"gene_exact:{gene.lower()} AND organism_id:{taxid} AND reviewed:true",
        f"gene:{gene.lower()}* AND organism_id:{taxid}",
        f"gene:{gene} AND organism_id:{taxid}",
        f"gene:{gene.lower()} AND organism_id:{taxid}",
    ]
    for query in queries:
        df = _get_uniprot_search_rows(query)
        if df.empty:
            continue
        preferred = df[df["Gene Names (primary)"].fillna("").str.lower().str.startswith(gene.lower())]
        hit = preferred.iloc[0] if not preferred.empty else df.iloc[0]
        status = "reviewed_search_hit" if str(hit.get("Reviewed", "")).strip().lower() == "reviewed" else "best_available_search_hit"
        return str(hit["Entry"]), status
    return None, "not_found"


def _download_fasta_for_accession(accession: str, dest: Path) -> SeqIO.SeqRecord:
    download(f"https://rest.uniprot.org/uniprotkb/{accession}.fasta", dest)
    return next(SeqIO.parse(dest, "fasta"))


def _alignment_maps(alignment: AlignIO.MultipleSeqAlignment) -> tuple[dict[str, str], dict[str, dict[int, int]]]:
    aligned = {record.id: str(record.seq) for record in alignment}
    seq_to_aln: dict[str, dict[int, int]] = {}
    for seq_id, seq in aligned.items():
        mapping = {}
        residue_idx = 0
        for aln_idx, aa in enumerate(seq):
            if aa != "-":
                residue_idx += 1
                mapping[residue_idx] = aln_idx
        seq_to_aln[seq_id] = mapping
    return aligned, seq_to_aln


def _conservation_status(column_aas: list[str], wt_aa: str) -> str:
    nongap = [aa for aa in column_aas if aa != "-"]
    if not nongap:
        return "unresolved"
    if all(aa == wt_aa for aa in nongap):
        return "identical"
    wt_group = AA_GROUPS.get(wt_aa, "other")
    if all((aa == wt_aa) or (AA_GROUPS.get(aa, "other") == wt_group) for aa in nongap):
        return "conservative"
    return "variable"


def _consurf_like_grade(paralog_fraction: float, ortholog_fraction: float) -> int:
    fractions = [x for x in [paralog_fraction, ortholog_fraction] if pd.notna(x)]
    if not fractions:
        return np.nan
    mean_fraction = float(np.mean(fractions))
    return int(max(1, min(9, round(mean_fraction * 9))))


def _clean_residue_label(row: pd.Series) -> str:
    return f"{row['aa_ref']}{int(row['residue_num'])}"


def compute_top10_paralog_conservation(df_top10: pd.DataFrame) -> pd.DataFrame:
    fasta_path = CACHE_DIR / "kcnq_paralogs.fasta"
    aln_path = CACHE_DIR / "kcnq_paralogs.aln.fasta"
    records = []
    for gene, uniprot in UNIPROT_IDS.items():
        dest = CACHE_DIR / f"{gene}_{uniprot}.fasta"
        download(f"https://rest.uniprot.org/uniprotkb/{uniprot}.fasta", dest)
        record = next(SeqIO.parse(dest, "fasta"))
        record.id = gene
        record.name = gene
        record.description = gene
        records.append(record)
    SeqIO.write(records, fasta_path, "fasta")
    subprocess.run(
        ["clustalo", "-i", str(fasta_path), "-o", str(aln_path), "--force", "--outfmt=fasta"],
        check=True,
        capture_output=True,
        text=True,
    )
    alignment = AlignIO.read(aln_path, "fasta")
    aligned, seq_to_aln = _alignment_maps(alignment)

    rows = []
    ortholog_accessions: dict[str, dict[str, str | None]] = {}
    ortholog_status: dict[str, dict[str, str]] = {}
    ortholog_alignments: dict[str, tuple[dict[str, str], dict[str, dict[int, int]]]] = {}
    for gene in UNIPROT_IDS:
        ortholog_accessions[gene] = {}
        ortholog_status[gene] = {}
        gene_records = []
        for species_label in ORTHOLOG_SPECIES:
            accession, status = _find_ortholog_accession(gene, species_label)
            ortholog_accessions[gene][species_label] = accession
            ortholog_status[gene][species_label] = status
            if accession is None:
                continue
            dest = CACHE_DIR / f"{gene}_{species_label}_{accession}.fasta"
            record = _download_fasta_for_accession(accession, dest)
            record.id = species_label
            record.name = species_label
            record.description = f"{species_label}|{accession}"
            gene_records.append(record)
        gene_fasta = CACHE_DIR / f"{gene}_orthologs.fasta"
        gene_aln = CACHE_DIR / f"{gene}_orthologs.aln.fasta"
        if gene_records:
            SeqIO.write(gene_records, gene_fasta, "fasta")
            subprocess.run(
                ["clustalo", "-i", str(gene_fasta), "-o", str(gene_aln), "--force", "--outfmt=fasta"],
                check=True,
                capture_output=True,
                text=True,
            )
            gene_alignment = AlignIO.read(gene_aln, "fasta")
            ortholog_alignments[gene] = _alignment_maps(gene_alignment)

    for row in df_top10.itertuples(index=False):
        residue_num = int(row.residue_num)
        aln_idx = seq_to_aln.get(row.gene, {}).get(residue_num)
        same_fraction = np.nan
        grade = np.nan
        aa_by_paralog = {gene: pd.NA for gene in UNIPROT_IDS}
        paralog_identity_n = pd.NA
        paralog_status = "unresolved"
        if aln_idx is not None:
            column = [aligned[gene][aln_idx] for gene in aligned]
            nongap = [aa for aa in column if aa != "-"]
            if nongap:
                wt_aa = str(row.aa_ref)
                same_fraction = sum(aa == wt_aa for aa in nongap) / len(nongap)
                grade = int(max(1, min(9, round(same_fraction * 9))))
                paralog_identity_n = int(sum(aa == wt_aa for aa in nongap))
                paralog_status = _conservation_status(column, wt_aa)
                for gene_name in aligned:
                    aa_by_paralog[gene_name] = aligned[gene_name][aln_idx]
        ortholog_fraction = np.nan
        ortholog_identity_n = pd.NA
        ortholog_status_label = "unresolved"
        ortholog_accession_summary = []
        ortholog_map = ortholog_alignments.get(row.gene)
        if ortholog_map is not None:
            orth_aligned, orth_seq_to_aln = ortholog_map
            orth_aln_idx = orth_seq_to_aln.get("human", {}).get(residue_num)
            if orth_aln_idx is not None:
                orth_column = [orth_aligned[species][orth_aln_idx] for species in orth_aligned]
                nongap_orth = [aa for aa in orth_column if aa != "-"]
                if nongap_orth:
                    wt_aa = str(row.aa_ref)
                    ortholog_identity_n = int(sum(aa == wt_aa for aa in nongap_orth))
                    ortholog_fraction = ortholog_identity_n / len(ORTHOLOG_SPECIES)
                    ortholog_status_label = _conservation_status(orth_column, wt_aa)
        for species_label, accession in ortholog_accessions.get(row.gene, {}).items():
            status = ortholog_status.get(row.gene, {}).get(species_label, "not_found")
            ortholog_accession_summary.append(f"{species_label}:{accession or 'NA'}({status})")
        grade = _consurf_like_grade(same_fraction, ortholog_fraction)
        rows.append(
            {
                "gene": row.gene,
                "protein_change": row.protein_change,
                "residue_num": residue_num,
                "candidate": f"{row.gene} {row.protein_change}",
                "residue_label": f"{row.aa_ref}{residue_num}",
                "aa_in_KCNQ1": aa_by_paralog["KCNQ1"],
                "aa_in_KCNQ2": aa_by_paralog["KCNQ2"],
                "aa_in_KCNQ3": aa_by_paralog["KCNQ3"],
                "aa_in_KCNQ4": aa_by_paralog["KCNQ4"],
                "aa_in_KCNQ5": aa_by_paralog["KCNQ5"],
                "paralog_identity_n": paralog_identity_n,
                "paralog_identity_label": f"{paralog_identity_n}/5" if pd.notna(paralog_identity_n) else pd.NA,
                "paralog_conservation_status": paralog_status,
                "paralog_conservation_fraction": same_fraction,
                "cross_species_identity_n": ortholog_identity_n,
                "cross_species_identity_label": f"{ortholog_identity_n}/4" if pd.notna(ortholog_identity_n) else pd.NA,
                "cross_species_conservation_fraction": ortholog_fraction,
                "cross_species_conservation_status": ortholog_status_label,
                "paralog_conservation_grade_1to9": int(max(1, min(9, round(same_fraction * 9)))) if pd.notna(same_fraction) else np.nan,
                "consurf_like_grade_1to9": grade,
                "ortholog_accessions": " | ".join(ortholog_accession_summary),
            }
        )
    return pd.DataFrame(rows)
