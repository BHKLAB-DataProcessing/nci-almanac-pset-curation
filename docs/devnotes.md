# Development decisions

See [processing decisions](implementation.md), [source provenance](data_sources.md) and [validation](validation.md). Dependencies and source checksums are pinned; generated data and diagnostic artifacts remain outside Git.

## NOTZ reproduction from public raw data

The production response is now `PERCENTGROWTHNOTZ / 100`, retaining that reported
value and using `pmax(0, viability_reported)` for observed calculations. No upper
cap is imposed on observed mono references. Curve fitting still truncates to
[0,1], matching the peer fitter's default behavior. Mono fitting and combination
summaries use sample + SCREENER_STUDY context, preserving NSC-specific identity.
Original physical observations remain in `raw`; derived scores and matched fits
are separate `synergy` and `profiles` assays. Stable source panel/cell IDs and
NSCs remain public keys, avoiding the peer object's sample-alias mismatch.

The peer script was built with PharmacoGx 3.7.5 at
3dab6c1e4f7d43b03398d80a4873f21e99bacef6. This pipeline uses the requested PR178
snapshot (3.14.1, 28c58707eece80086bee75d8ef0e8afae96389ca), explicitly setting
`viability_as_pct=FALSE` for fit-derived metrics. The peer omitted that flag for
IC50 and its historically named computeAUC call, so its fitted values are not
assumed numerically equivalent. AUC and AAC remain distinct here; observed HSA
and Bliss and raw-source values are the primary exact-concordance targets.

The saved peer PSet is held out for validation, never read by production stages.
All 3,686,475 peer source observations match the public November 2017 CSV across
the checked identifiers and raw measurements. Molecular assays are rebuilt from
the documented public CellMiner/DTP inputs, not copied from NCI60_2021.qs.
The old 74 fit checkpoints and pre-NOTZ response are preserved. Changed input
and settings hashes create a separate fit cache. The output is versioned as
`NCI_ALMANAC_PSet.rds`. Full source-observation and final-object validation passed; fitted-reference differences remain documented in `validation.md`.

### Observed-well indexing

CoreGx 2.7.3 builds a Cartesian row/column index before removing unobserved cells, exceeding the 16 GiB R vector limit for the full Matrix screen. The repo-local `construct_observed_tre.R` adapter constructs the equivalent raw-assay index from observed physical wells. It preserves the CoreGx object layout and is checked against the standard constructor, including missing identifiers, repeated doses, and subsequent assay assignment. No response measurements or fit definitions change. Recheck equivalence when upgrading CoreGx.

`pixi run setup` now runs `workflow/tests/test_observed_tre.R` after package installation. This checks exact object equivalence to the installed CoreGx public constructor and assay assignment, so a future dependency change cannot silently bypass this compatibility check.

### Molecular assay type annotations

Assembly sets each molecular experiment's `metadata$annotation` to its explicit assay type; each assay retains its correct molecular type. This is required by the current MCP `drugSensitivitySig` consumer. Source values and units are unchanged. Proteomics, methylation and heterogeneous DTP experiments retain their true types even though that consumer does not currently support all of them. Use explicit numeric `sProfiles` for combination-response analyses; do not delete descriptive profile fields merely to make a consumer run.

Package-description installation paths (`file` attributes and `LibPath`) are removed from the public object annotation. Package versions and GitHub revisions remain intact, and full local session/package manifests are retained outside the RDS under `data/procdata/metadata`. QC checks that the serialized annotation is already portable.

## Full NOTZ source-observation validation

An independent held-out check reconstructed mono references from the public-source response import using sample, SCREENER/STUDY and drug/dose context. All 2,873,514 combination observations matched Nikta’s reference object exactly for NOTZ modeling viability, both observed mono references, HSA and Bliss, including missing-value patterns. Maximum absolute difference was zero for all five quantities. This comparison did not use the peer object as a production input and did not run or validate the new fitted metrics. Evidence: `data/procdata/metadata/peer_notz_full_source_observation_concordance.tsv`. Full fitted output, assembly and reloaded QC remain separately required.

## Held-out comparison criteria

Observed NOTZ viability, observed mono references and observed HSA/Bliss require matching context keys, equal missingness patterns and maximum absolute error below 1e-10. The independent full-source check already observed zero difference for all 2,873,514 combination contexts. Final exported tables are checked again after assembly.

Fitted metrics are compared separately at both combination-profile drug-1 level and deduplicated monotherapy context level combining both drug roles. The latter prevents repeated combination profiles from over-weighting a mono fit. Conflicting peer fits within one NSC/sample/technical context are rejected. Diagnostics report finite coverage, missingness differences, median/p90/p99/max absolute differences and rank correlation. Threshold counts are descriptive, not equivalence gates: 1e-8 on each reported scale, a twofold IC50 difference (absolute log10 difference <= log10(2)), and AAC absolute difference <= 0.01. No undocumented fit-equivalence threshold is applied.

The peer's historical computeAUC field uses complementary response-area semantics on a 0–100 scale; compare it divided by 100 with current computeAAC, not current viability-area computeAUC. The peer's omitted viability_as_pct arguments and older Hill model can cause substantive fitted differences. These cannot be dismissed as roundoff or implied equivalent by the observed-metric match.
