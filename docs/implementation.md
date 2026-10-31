# NCI ALMANAC processing decisions

The public ALMANAC growth table supplies percent growth. Source percent growth is preserved; viability for the reference workflow is max(0, percent growth / 100). Single-agent curves are pooled within screener/study, including matching wells on different plates. CellMiner and DTP supply publicly processed molecular measurements, following the NCI60 pipeline; no peer NCI60 PSet is required.

## Package and measurement conventions

- `packages.PharmacoGx` selects [PR #178](https://github.com/bhklab/PharmacoGx/pull/178), branch `feat/mcp`. Validation installed commit `28c58707eece80086bee75d8ef0e8afae96389ca` (3.14.1); CoreGx main resolved to `bf6d814abd3dc0fc98f9f419f430cc6473046cb2` (2.7.3).
- `renv.lock` pins the tested R dependencies, Bioconductor release and exact CoreGx/PharmacoGx commits. `pixi.lock` pins the system and R environment; setup restores into the project’s Pixi library and verifies the resolved package versions.
- PR178 Hill fits retain `E0`, `E_inf`, slope, EC50, IC50, fit status and Rsquare. Fits use responses clipped to [0,1]; source measurements remain available. Doses retain their source unit in the TRE.
- AUC means normalized **viability** area under PR178; AAC means normalized response area. Both are exported under distinct names.
- HSA/Bliss observed and recomputed fitted references are separate. Fitted predictions use the fitted E0. Positive excess is reference viability minus combination viability. Missing monotherapy references remain missing.
- The package's Loewe and ZIP APIs assume E0=1. Comparison columns explicitly use the `Loewe_E0_1` / `ZIP_E0_1` prefix and must not be interpreted as predictions from the complete four-parameter fit. They are computed for eligible Matrix/ALMANAC pairs, following the peer API path. Negative/nonfinite slopes or invalid doses/EC50 are excluded.
- Source IDs remain primary keys. AnnotationDB enriches cell lines and compounds; ambiguous matches remain unmapped. AnnotationGx fields are NA because that method was not run. Gene mapping uses the configured AnnotationDB-hosted GENCODE v50 reference; gene `In.AnnotationDB` denotes membership in that reference, not a separate API lookup.
- Unique current Ensembl IDs are used for gene-level assays. Multiple features mapping to a gene use lexical-first source feature, with an exclusion audit. Fusions retain paired gene IDs. Heterogeneous DTP targets retain source feature IDs and separate method/unit assays.
- Molecular `batchid` identifies the source processing cohort, not an inferred experimental batch. Variant zeros mean no variant reported in the source, not confirmed wild type.

## Implementation and validation

The scripts implement download, extraction, response processing, annotation, curve fitting, synergy, molecular import, PSet construction and QC. A script being implemented does not imply that its full dataset run passed. See `validation.md` for the latest run state. Generated logs, fit checkpoints, receipts and mapping audits live under `data/`, outside Git.

Downloads are atomic, checksum-verified and resumable. Original archives are retained. Extraction rejects traversal/symlink entries and streams only gene records from GTF to avoid a multi-gigabyte uncompressed annotation. Expensive fits checkpoint complete chunks. A failed fit is retained with its error status; all-failed runs stop.

Four repositories are retained because response schemas and normalization differ substantially. Small utilities remain repo-local so each repository can run independently. No peer-created RDS or prebuilt PSet is a production input. PSet data and rendered reports are distributed separately from this source repository; no ORCESTRA upload has been performed.

## Executable source cleanup

Reference-only Anchor portal/fitted/day-1 tables and ALMANAC field documentation were removed from the production `sources` list. Their existing downloads and source links remain in the historical provenance notes. Matrix portal tables remain required for compound/model metadata.

## Four-parameter dose support

Three distinct doses preserve the peer fitting threshold, but cannot independently identify four Hill parameters. QC reports `four_parameter_dose_support` and dose counts; successful optimization is not evidence that the fit is identifiable. ALMANAC observed monotherapy means are floored at zero without an upper cap; fitting clips to [0,1], while raw percent growth remains intact.
