# NCI ALMANAC PharmacoSet curation

## Overview

Minimal ADAPT PSet curation template. Sources are pending review; processing is intentionally unimplemented.

## Quick Start

1. Review `config/pipeline.yaml` and [source notes](docs/data_sources.md).
2. Select releases, direct-download URLs, checksums, metadata and molecular sources.
3. When ready for implementation, use `pixi install --locked`, then `pixi run setup`.

Setup installs the latest GitHub default-branch HEAD of CoreGx and PharmacoGx explicitly into the Pixi R library. It records resolved SHAs in `data/procdata/metadata/github_packages.tsv`. Pin those refs before a reproducible production run. GitHub R packages and their transitive R dependencies are outside the Pixi lockfile; setup has not been run.

## Pipeline

Planned: configurable raw downloads → processing/normalization → harmonization → PharmacoSet → RDS and secondary TSV tables.

`pixi run pipeline` invokes `workflow/Snakefile`; its placeholder curation script deliberately stops without writing a PSet. No download rules or dataset processing have been implemented. Expected end-to-end runtime: unknown until sources and fitting are selected.

## Inputs

Source candidates live in `config/pipeline.yaml`. Null URLs/releases are unresolved decisions, not working inputs. Legacy scripts are references only; prebuilt PSets will not substitute for raw data.

## Outputs

- `data/results/NCI_ALMANAC_PSet.rds` — canonical PharmacoSet, written with `saveRDS()` after implementation.
- `data/results/tables/*.tsv` — planned secondary access tables; `pixi run zipup` archives them with maximum gzip compression.
- `data/procdata/metadata/` — future generated annotations, source/software manifests and mapping misses.

RDS only for the serialized PharmacoSet. Generated data are ignored by Git.

## QC and Docs

`workflow/qc/pset_qc.Rmd` is a report outline. `pixi run qc` will render to `docs/pset_qc.html` once implemented; it currently stops explicitly. See [sources](docs/data_sources.md) and [decisions](docs/devnotes.md).

## Development

Maintain dependencies with `pixi add`; commit `pixi.toml` and `pixi.lock`. Keep dataset parameters in `config/pipeline.yaml`, logic in `workflow/scripts/`, source metadata in `data/rawdata/metadata/`, and generated metadata in `data/procdata/metadata/`. Split a larger DAG under `workflow/rules/` when needed. Do not add a separate logs directory. CI checks manifest syntax and tracked-data hygiene only; no pipeline runs.
