# Contributing

This repository is organized as a publication companion rather than a general-purpose software package, so contributions should preserve reproducibility of the frozen release.

## Guidelines

- Keep code changes modular inside [kcnq_pipeline](kcnq_pipeline) or [scripts](scripts).
- Do not overwrite the frozen outputs in [results_final](results_final) without documenting the rerun conditions and regenerated files.
- Keep large downloaded resources, caches, and temporary docking intermediates out of version control.
- When benchmark definitions or thresholds change, update the corresponding export tables and the repository README together.

## Reporting Issues

Please include:

- the script or module involved
- the exact command used
- the relevant input file paths
- the traceback or failing output

## Pull Requests

Pull requests should explain:

- what changed
- whether the frozen outputs were regenerated
- which figures or export tables are affected
