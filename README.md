# KCNQ Rescue Prioritization Pipeline

Computational pipeline for family-wide prioritization of potentially rescuable missense variants in `KCNQ1-5`, with structural annotation, composite pathogenicity scoring, tractability-aware ranking, wild-type docking, and benchmarking against published functional datasets.

This repository is the cleaned, publication-facing codebase derived from the exploratory notebook workflow archived under [archive/legacy_notebook_run](archive/legacy_notebook_run). The reproducible pipeline lives in [kcnq_pipeline](kcnq_pipeline) and [scripts](scripts); frozen publication outputs live in [results_final](results_final).

## Highlights

- Integrates ClinVar, gnomAD, CADD, and AlphaMissense into a composite pathogenicity layer.
- Annotates missense variants across `KCNQ1-5` using curated structures for `KCNQ1`, `KCNQ2`, `KCNQ4`, `KCNQ5`, and AlphaFold for `KCNQ3`.
- Scores each clinical variant with a structural opportunity score and a tractability-adjusted rescue priority score.
- Runs wild-type docking for prioritized candidates and cross-drug docking matrices with Vina.
- Benchmarks against Bhatt `KCNQ2`, Vanoye `KCNQ1`, and Brewer `KCNQ1` functional datasets.
- Exports publication-ready figures, tables, and a consolidated workbook.

## Repository Layout

- [kcnq_pipeline](kcnq_pipeline): pipeline modules for fetching, structural annotation, scoring, docking, benchmarking, and figure generation.
- [scripts](scripts): staged entry points used to reproduce the frozen release.
- [pipeline_data/benchmarks](pipeline_data/benchmarks): benchmark spreadsheets and processed benchmark tables used in the manuscript analyses.
- [pipeline_data/structures](pipeline_data/structures): receptor structure files used for annotation and docking.
- [results_final/exports](results_final/exports): frozen analysis tables, benchmark summaries, and final workbook.
- [results_final/figures](results_final/figures): frozen main and supplementary figures.
- [archive](archive): provenance material from the earlier notebook-driven workflow and internal drafting artifacts.

## Environment

The pipeline was prepared against Python `3.11`.

Install the Python dependencies:

```bash
python -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
```

Additional tooling is required for docking:

- `RDKit`
- `AutoDock Vina`
- `Open Babel` (`obabel`)
- `mk_prepare_ligand.py` from Meeko/AutoDockTools

The docking stage is the only part of the pipeline that requires these external tools. The scoring, benchmarking, export, and figure stages can run from the frozen outputs without rerunning docking.

## Reproducing the Frozen Analysis

Run the staged scripts in this order from the repository root:

```bash
python scripts/01_fetch_variants.py
python scripts/02_structural_annotation.py
python scripts/04_scoring.py
python scripts/03_docking.py
python scripts/05_benchmarking.py
python scripts/06_figures.py
python scripts/07_export.py
```

Important notes:

- [pipeline_data/raw](pipeline_data/raw) is intentionally ignored because it can contain large downloaded resources such as the AlphaMissense genome-wide file.
- [scripts/fetch_revel_baseline.py](scripts/fetch_revel_baseline.py) is the only optional network-dependent benchmark helper not run by default in the staged pipeline.
- The repository keeps the frozen final outputs under [results_final](results_final) so the published figures and tables remain reproducible without recomputing every network call.

## Main Outputs

Primary release artifacts:

- [kcnq_rescue_analysis_final.xlsx](results_final/exports/kcnq_rescue_analysis_final.xlsx)
- [kcnq_top10_rescue_candidates_final.csv](results_final/exports/kcnq_top10_rescue_candidates_final.csv)
- [kcnq_all_annotated_variants_final.csv](results_final/exports/kcnq_all_annotated_variants_final.csv)
- [fig1_umap.png](results_final/figures/fig1_umap.png)
- [fig2_structural_annotation.png](results_final/figures/fig2_structural_annotation.png)
- [fig3_wt_docking.png](results_final/figures/fig3_wt_docking.png)
- [fig4_docking_matrix.png](results_final/figures/fig4_docking_matrix.png)
- [fig5_top10_ranked.png](results_final/figures/fig5_top10_ranked.png)
- [fig6_benchmarking.png](results_final/figures/fig6_benchmarking.png)
- [fig7_pose_renders.png](results_final/figures/fig7_pose_renders.png)

## Current Benchmark Summary

These values are drawn from the frozen tables in [results_final/exports](results_final/exports):

- `KCNQ2` Bhatt retigabine benchmark, structural opportunity score:
  - ROC-AUC `0.921` for `responsive_any`
  - ROC-AUC `0.889` for `responsive_current50`
  - ROC-AUC `0.283` for `reaches_wt_current`
- `KCNQ2` Bhatt retigabine benchmark, rescue priority score:
  - ROC-AUC `0.873` for `responsive_any`
  - ROC-AUC `0.825` for `responsive_current50`
  - ROC-AUC `0.654` for `reaches_wt_current`
- `KCNQ1` Vanoye overlap:
  - Spearman `rho = -0.374`
  - `p = 0.00411`
- Structural validation across the clinical set:
  - Mann-Whitney `p = 2.15e-80`
  - Cohen's `d = -1.07`
  - logistic `OR per Å = 0.9555`

For the complete benchmark and ablation inventory, see the CSV exports in [results_final/exports](results_final/exports).

## Datasets 

Benchmark spreadsheets bundled with the repository:

- Bhatt `KCNQ2` supplemental tables
- Vanoye `KCNQ1` dataset
- Brewer `KCNQ1` MAVE tables

Large external downloads are not committed. The largest omitted file is the genome-wide AlphaMissense TSV, which should be placed under [pipeline_data/raw](pipeline_data/raw) if a full re-run is required.

## Citation

If you use this code or the frozen outputs, cite the repository metadata in [CITATION.cff](CITATION.cff). If a manuscript DOI becomes available, update the `preferred-citation` block there.

## License

Released under the [MIT License](LICENSE).
