# NCI ALMANAC data sources

Status: source review pending. This is a candidate registry, not an approved input manifest. `config/pipeline.yaml` is the source of truth for URLs and selected versions.

| Source | Version/release | URL / access | Format |
|---|---|---|---|
| NCI ALMANAC combination growth | November 2017 | https://wiki.nci.nih.gov/spaces/NCIDTPdata/pages/338237347/NCI-ALMANAC | CSV |
| NCI60 raw molecular data and metadata | Unresolved | Unresolved | To be selected |

Direct Figshare URLs, where present, were recovered from historical local article metadata and have not been revalidated in this task. Confirm availability, licensing/access, exact release, raw-versus-derived status, expected files and checksums before adding download rules. HMCL and NCI direct URLs remain unset until exact assets are selected. Automatic direct downloads should be the first DAG stage wherever public access permits.

Review growth/viability semantics, technical replicate grouping and synergy recomputation. The legacy script imports PSet_NCI60.rds for molecular profiles; identify downloadable raw molecular/metadata sources instead of carrying that prebuilt PSet into this raw-data pipeline.

Expected generated output: `data/results/NCI_ALMANAC_PSet.rds` plus secondary TSV exports. No input files or generated PSets have been copied into this repo.

Legacy processing references (inspect before porting):
- `/Users/michael/Projects/BHKLab/PGx/pgx3/Nasim/case_study3_nci_almanac_curation (5).R`
