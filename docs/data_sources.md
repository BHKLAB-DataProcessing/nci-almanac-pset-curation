# NCI ALMANAC data sources

Response measurements come from the original [NCI ALMANAC](https://wiki.nci.nih.gov/spaces/NCIDTPdata/pages/338237347/NCI-ALMANAC) November 2017 growth export: 3,686,475 rows for 60 lines and 105 NSCs. Both PERCENTGROWTH and PERCENTGROWTHNOTZ are retained. Reported viability is PERCENTGROWTHNOTZ/100; observed calculations floor at zero and fitting additionally caps at one. Molecular profiles are imported directly from the configured public CellMiner and NCI DTP processed datasets. No previous NCI60 PSet is required. Molecular releases are not asserted identical to historical reference PSets.

## Source files

| File | Role | Bytes | Download URL |
|---|---|---:|---|
| `ComboDrugGrowth_Nov2017.zip` | published screening archive | 85966226 | [source](https://wiki.nci.nih.gov/download/attachments/338237347/ComboDrugGrowth_Nov2017.zip?version=1&modificationDate=1510057275000&api=v2) |
| `ALMANAC_DataFields.txt` | source field and compound reference | 3880 | [source](https://wiki.nci.nih.gov/download/attachments/338237347/ALMANAC_DataFields.txt?version=1&modificationDate=1513947309000&api=v2) |
| `ComboCompoundNames_all.txt` | source field and compound reference | 50857 | [source](https://wiki.nci.nih.gov/download/attachments/338237347/ComboCompoundNames_all.txt?version=1&modificationDate=1493822512000&api=v2) |

## Executable sources and annotation

`config/pipeline.yaml` is the authoritative list of production/QC downloads and expected checksums. Reference-only files in the table need not be downloaded by production stages. Generated source receipts, mapping exclusions and API responses are stored under ignored `data/` directories. AnnotationDB live mode can change independently of measured data; replay mode verifies reviewed response hashes and fails on missing or changed responses. See the README for the exact replay limitation. Gene identifiers use the configured GENCODE reference; mapping does not realign or renormalize source measurements.

## Molecular release and scales

Each source workbook records CellMiner database version **2.15**, human genome **HG-19**, and the source date string
**`09-11-2025`** (preserved verbatim rather than interpreting the date order).
These are upstream source metadata; the PSet build date is not the source release
date. Gene identifier mapping to the configured current annotation does not
realign source reads or change their original genome assembly.

| Input ID | Source value scale | Pipeline treatment |
|---|---|---|
| `rnaseq_composite` | `log2(FPKM + 1)` | Retain downloaded values; no second log transform |
| `cnv_gene_summary` | Average log2 intensity ratio, cell-line DNA versus normal DNA | Retain ratio; do not label as absolute copy number |
| `methylation_gene_average` | Beta values, 0 unmethylated to 1 methylated | Retain gene-average beta values |
| `exome_protein_function_affecting` | Summation of variants | Retain counts; not automatically binary mutation status |
| `swath_protein` | Source SWATH values | Retain source abundance; do not infer a new transformation |


## Reproducible naming

The main curation derives NSC/formulation names and CAS annotations from the original NCI `ComboCompoundNames_small.txt` and `ComboCompoundSet.sdf`, whose URLs and checksums are configured. `config/naming-policy.json` contains reviewed display-label exceptions, not response or molecular data. Sample names originate in the screen. NSC 122758 is Tretinoin; NSC 125066 remains the Bleomycin sulfate product; NSC 721517 is Zoledronic acid. Vemurafenib NSCs 753082 and 761431 share a display label, with both source identifiers retained in response keys and treatment metadata. No chemical equivalence or pooling of measurements is implied.

Naming is part of `curate.R` and requires no historical PSet or local NCI-60 cache. RNA-seq is exposed as `rnaseq.comp`/`exprs`; values are unchanged. The explicit context-level biomarker accessor is in `workflow/scripts/ALMANAC_biomarker_helpers.R`.
