# NCI ALMANAC PharmacoSet curation

## Overview

An ADAPT curation pipeline from public downloadable sources to a PharmacoSet RDS. The full current artifact and source-to-object QC checks passed; [validation status](docs/validation.md) distinguishes complete checks from pending work. [Processing decisions](docs/implementation.md) explain the reference-pipeline correspondence and source boundaries.

## Quick Start

```sh
pixi install --locked
pixi run setup
pixi run pipeline
pixi run qc
pixi run zipup
```

Review `config/pipeline.yaml` first. `packages` records fixed GitHub revisions; `sources` records original download URLs, IDs and checksums. Setup installs into the dedicated Pixi R library and records resolved revisions. It refuses a global R library. Contact: michaelcao-anh.tran@uhn.ca.

## Pipeline

`workflow/Snakefile` orchestrates source downloads → safe extraction → response/molecular processing → annotation → fitting and synergy → PharmacoSet construction and TSV exports. `pixi run qc` separately renders the report. Fit checkpoints allow restart. Default fitting uses four workers with one concurrent pipeline job; reduce workers if memory is constrained.

Full end-to-end runtime has not yet been measured. Downloads, curve fitting and molecular imports dominate. Allow several GiB for source archives, intermediate objects and results.

## Inputs

The public ALMANAC growth table supplies both PERCENTGROWTH and PERCENTGROWTHNOTZ. Both source fields are preserved; this version follows Nikta’s latest script with reported viability = PERCENTGROWTHNOTZ/100 and modeling viability = max(0, reported viability). The peer PSet is held out for comparison only and is never a production input. Single-agent curves are pooled within screener/study, including matching wells on different plates. CellMiner and DTP supply publicly processed molecular measurements, following the NCI60 pipeline; no peer NCI60 PSet is required.
Full response import checked: 3,686,475 dose rows; 60 stable panel/cell IDs; 105 NSC compounds. Source counts are not annotation match counts. Original inspection records and links remain in [data sources](docs/data_sources.md); `config/pipeline.yaml` is the executable source list. SHA-256 receipts bind a run to actual downloaded bytes, including mutable URLs.

## Outputs

- `data/results/NCI_ALMANAC_PSet.rds`: validated serialized PharmacoSet.
- `data/results/tables/*.tsv`: sample, treatment, sensitivity, curation, fit, combination and molecular-summary tables.
- `data/procdata/metadata/`: package revisions, session info, source receipts, annotation misses, feature exclusions and object validation.
- `docs/pset_qc.html`: rendered QC after a successful pipeline run.

The canonical PSet uses `saveRDS()`. Data and generated reports are ignored by Git. `zipup` creates a maximum-gzip TSV archive.

## QC and Docs

The single full-data QC report explores response distributions, fit support, representative curves, molecular coverage and small dataset-appropriate analyses. Published comparisons use matched metrics, contexts and units, with missingness, Pearson/Spearman correlations, agreement/bias and outlier diagnostics. See [QC standard](docs/qc-standard.md), [implementation](docs/implementation.md), [validation](docs/validation.md), and [source history](docs/data_sources.md). Lightweight integration fixtures remain internal software tests.

## Development

Maintain Conda dependencies with `pixi add` and retain `pixi.lock`. R dependencies installed from GitHub/Bioconductor are recorded separately. Keep processing in `workflow/scripts/`, QC in `workflow/qc/`, and generated files under `data/{rawdata,procdata,results}`. Never replace public input provenance with an opaque peer PSet.


### Canonical names and source identities

The main workflow derives complete maps from the original NCI compound-name/structure files and source sample metadata, then applies the reviewed display labels in `config/naming-policy.json`. It writes `data/results/NCI_ALMANAC_PSet.rds`: 104 named treatments retain all 105 NSCs, including separate Vemurafenib source contexts. Response values and fits are unchanged by naming. The RNA-seq experiment is `rnaseq.comp` with inner assay `exprs`. The explicit biomarker helper remains available for the context-level profiles.

## Rebuilding from original sources

Run `pixi install --locked`, `pixi run setup`, then `pixi run pipeline`. `pixi.lock` fixes the system/R environment; `renv.lock` fixes the tested R dependency versions, Bioconductor release and CoreGx/PharmacoGx commits. Setup modifies only this project’s Pixi library. Input data, API responses, fitted intermediates and final RDS/HTML files are generated and ignored by Git. No existing PSet is required.

Original scientific downloads are checksum-pinned in `config/pipeline.yaml`; a changed remote file stops the build. AnnotationDB defaults to `mode: live`, reusing local responses when available. Set `refresh: true` for an explicit live refresh. API annotations can change independently of the measured data. `mode: replay` with `refresh: false` requires every response to match `config/annotationdb-checksums.json`; it never silently accepts changed bytes. Those checksums describe the reviewed cache, not a publicly archived API snapshot. A fresh replay can fail if the service changes, in which case the original response cache is needed; a live build remains possible. Generated response hashes are recorded in `data/procdata/metadata/annotationdb_responses.tsv`.

Verification covered a clean macOS ARM64 installation of the identical shared dependency lock, each dataset’s existing setup checks, all four isolated workflow graphs, access to every configured source URL, and selected fresh download/extraction stages. The fresh HMCL response stage reproduced all 904,376 observations and its sample, treatment and published-profile tables exactly. Full fitting and end-to-end regeneration of all four PSets were not repeated; current final objects retain the previously verified fitted measurements.
