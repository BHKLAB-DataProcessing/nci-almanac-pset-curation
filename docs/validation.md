# Validation status

The canonical named PharmacoSet is rebuilt from the source-curated fitted intermediates. Its scientific fingerprints, sample/treatment identities, source NSCs and molecular experiments match the preceding named release exactly. The source-driven workflow now generates the naming maps directly from NCI downloads and reviewed label rules. No refit or potency patch is included. Publication requires acknowledging the fitted-metric limitations below.

## Artifact and completed checks

- `data/results/NCI_ALMANAC_PSet.rds`: 438,499,014 bytes; SHA-256 `ccc32392a525b8dcdc7fb5f8a800c743db82dc2aecc65b8d15babcbccacfc587`.
- 60 response samples, 104 named treatments retaining all 105 NSCs, 3,686,475 physical measurements, 2,873,514 derived combination observations, 260,263 monotherapy contexts and 311,604 combination profiles; 104 molecular assays.
- All 25 source-to-object checks pass: exact row counts, unique dose/replicate/well keys, identifiers, retained source measurements and missing values. This explicitly covers the CoreGx entity-only constructor-key collapse hazard; S4 validity alone was not used as evidence of preservation.
- Independent reload passes object validity, foreign keys, physical-well uniqueness, observed HSA/Bliss formulas and indexed-subset API checks. Naming retains distinct source NSCs in response and fit-context keys.
- All 2,873,514 peer combination observations align by NSC, sample, technical context and dose. Observed viability, both mono references, HSA and Bliss agree within floating-point roundoff (maximum absolute difference 5.56e-16), including matching missingness.
- The peer exposes 193,369 unique mono contexts through its combination profiles. Of these, 193,323 map to current source-based monotherapy contexts. The remaining 46 contain only missing fit metrics and have zero source monotherapy wells; they remain explicitly audited missing placeholders, not fabricated fits. Another 66,940 current source mono contexts are not exposed in those peer profiles.

Evidence is in `data/procdata/metadata/source_roundtrip.tsv`, `peer_notz_concordance.tsv`, `peer_monofit_context_coverage.tsv`, `peer_unmatched_empty_monofit_contexts.tsv` and the validation logs. The substantive report is [pset_qc.html](pset_qc.html), with full tables and figures in `docs/qc-data`.

## Scientific findings and limits

Observed-metric agreement does **not** establish fitted equivalence. Across unique monotherapy contexts:

| Metric | Finite paired contexts | Pearson | Spearman | Mean absolute difference | Other agreement summary |
|---|---:|---:|---:|---:|---|
| log10 IC50 (M) | 152,447 | 0.3872 | 0.6755 | 1.2510 log10 units | 29.70% within twofold; current-minus-peer mean bias -1.1708 log10 units |
| AAC fraction | 186,161 | 0.0987 | 0.0714 | 0.1802 | 13.50% within 0.01; current-minus-peer mean bias -0.1434 |

The peer used PharmacoGx 3.7.5 at `3dab6c1e4f7d43b03398d80a4873f21e99bacef6`, including historical complementary-area `computeAUC` on a 0–100 scale. Its saved area is divided by 100 for comparison with current AAC. The old script omitted `viability_as_pct=FALSE` for derived metrics; the current pipeline sets it explicitly and uses the requested MCP model. These are meaningful methodological differences. The observed poor agreement is not dismissed as numerical noise, and the report does not claim the current fitted summaries reproduce the peer's values.

Of 260,263 mono contexts, 247,145 successful fits use three usable doses, 5,936 use six, and 7,182 have insufficient usable doses. Three concentrations cannot independently identify four free Hill parameters. High in-sample R-squared is therefore not potency validation. The tested `computeIC50` default is relative half-effect concentration, not necessarily 50% absolute viability.

The pinned fitter also computes its reported variance-based R-squared before returning clipped asymptotes. In 69,628 contexts that statistic differs from the returned-curve diagnostic; 9,153 additional contexts fall below 0.5 after reevaluating the returned curve. The report preserves the original statistic and separately reports returned-curve variance-based and bias-sensitive SSE-based diagnostics. Extreme negative historical R-squared values are retained in statistics/outlier tables, with an explicitly labeled display zoom for readability.

A small matched-dose reproducibility analysis gives AAC Pearson 0.9191 and MAE 0.02834 across 390 finite pairs of distinct SCREENER/STUDY contexts. This is descriptive technical-context agreement, not independent biological validation. Observed-versus-fitted mean Bliss has Pearson 0.7257 and mean bias -0.03068 over 300,091 finite contexts; the choice of mono reference matters.

AnnotationDB maps 52/60 samples and 88/104 named treatments (89/105 source NSCs). Unmapped identifiers remain valid source keys and are exported for annotation review. The five CellMiner molecular layers cover 59 response samples; MDA-MB-468 is missing. The source label MDA-N is not silently equated with it. Heterogeneous DTP measurements remain 99 separate assays; family summaries describe coverage only.

## Reproducibility and release checks

The peer RDS is held out for validation and is never a production input. Original NCI measurements and molecular inputs come from configured public downloads. NOTZ response and checkpoint hashes are separate from the old 74 PERCENTGROWTH checkpoints, which remain preserved and incompatible. Molecular layers come from documented public CellMiner/DTP releases, not the historical peer NCI60 object.

`pixi.lock` and `renv.lock` pin R 4.5.3, Bioconductor 3.22 and all tested package versions. PharmacoGx 3.14.1 is pinned to `28c58707eece80086bee75d8ef0e8afae96389ca`; CoreGx 2.7.3 is pinned to `bf6d814abd3dc0fc98f9f419f430cc6473046cb2`. Setup is restricted to the project Pixi library. Scientific downloads are checksum-pinned; live AnnotationDB refresh and exact response replay are explicit separate modes. See the README for replay's archival limitation.

Clean exported code has passed dependency-graph validation and bounded original-source download/extraction checks. The shared locked environment was installed independently on macOS ARM64 and all four dataset setup checks passed. HMCL response preparation from a fresh source download reproduced all source-derived tables exactly. This does not constitute a fresh end-to-end rebuild of all four PSets: the current ALMANAC assembly reused verified fitted intermediates. No fitting or publication was performed during this regeneration repair.
