# Development decisions

- Decision: scaffold only; review sources before implementing downloads, normalization, fitting or serialization.
- Decision: separate HMCL, NCI ALMANAC, GDSC Matrix and GDSC Anchor repos. Matrix/Anchor share setup and mapping concepts but have different source schemas and processing. Keep helpers repo-local initially; revisit consolidation once implementation reveals stable overlap.
- Decision: canonical output is a portable PharmacoSet RDS, with optional TSV exports. Do not embed repository paths, download configs or source documents in public object annotation slots.
- Decision: install current GitHub CoreGx/PharmacoGx using explicit setup, record resolved SHAs and pin before production. No legacy feature branches or package-description patches are inherited. No package installation has been attempted.
- Planned harmonization: stable source primary keys; Cellosaurus sample mappings, PubChem compound mappings, and current/pinned GENCODE Ensembl gene IDs when molecular profiles are included. Keep biological attributes in TRE rowData/colData, not repeated in physical assay tables.
- Planned audits: export compound/cell/gene AnnotationDB mapping misses as TSV with source keys, relevant accessions, `In.AnnotationDB`, `AnnotationGx.Mapped` and `AnnotationGx.Title`. Keep fitting choices and API hosts in versioned config.
- Scaffold scope: no Copier migration, documentation publishing, release automation, containers, full QC, raw downloads or dataset implementation yet. These are deferred intentionally to keep this source-review template minimal.

## Dataset-specific review

Review growth/viability semantics, technical replicate grouping and synergy recomputation. The legacy script imports PSet_NCI60.rds for molecular profiles; identify downloadable raw molecular/metadata sources instead of carrying that prebuilt PSet into this raw-data pipeline.
