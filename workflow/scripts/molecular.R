# Source readers adapted from the NCI60 and GDSC public curation pipelines.
# No prebuilt peer PharmacoSet is an input.
source("workflow/scripts/common.R")
suppressPackageStartupMessages({library(SummarizedExperiment); library(MultiAssayExperiment)})
source("workflow/scripts/gene_helpers.R")
samples <- readRDS("data/procdata/response.rds")$sample
profiles <- list()

if (cfg$dataset$id == "HMCL") {
  x <- fread(source_file("screen", "/mutations/hmcl_predicted_mutations.tsv$"))
  x <- x[cell_line %in% samples$sampleid & !is.na(SNPEFF_GENE_NAME)]
  wide <- dcast(x, SNPEFF_GENE_NAME ~ cell_line, fun.aggregate=function(x) as.integer(length(x)>0), value.var="SNPEFF_GENE_NAME", fill=0L)
  matrix <- as.matrix(wide[, -1]); rownames(matrix) <- wide[[1]]
  profiles[["mutation.exome"]] <- gene_se(matrix, wide[[1]], wide[[1]], "mutation.exome", "binary source-reported variant presence")
}

if (cfg$dataset$id == "NCI_ALMANAC") {
  sample_lookup <- data.table(sampleid=samples$sampleid, sample_key=normal_key(samples$sample_name))
  sample_lookup[, sample_key := sub("ATCC$", "", sub("TB$", "", sample_key))]
  if (anyDuplicated(sample_lookup$sample_key)) stop("Ambiguous normalized NCI sample names")
  assay_types <- c(rnaseq_composite="rnaseq.composite", cnv_gene_summary="cnv.acgh",
                   methylation_gene_average="methylation.450k", exome_protein_function_affecting="mutation.exome",
                   swath_protein="proteomics.swath")
  units <- c(rnaseq_composite="log2(FPKM + 1)", cnv_gene_summary="log2 copy-number ratio",
             methylation_gene_average="gene-average beta", exome_protein_function_affecting="protein-function-affecting variant count",
             swath_protein="source SWATH abundance")
  for (id in names(assay_types)) {
    raw <- as.data.table(readxl::read_xls(source_file(id, "/output/[^/]+\\.xls$"), col_names=FALSE, col_types="text", .name_repair="minimal"))
    header <- which(grepl("^(Gene name|Identifier)(\\s|$)", raw[[1]]))[1]
    if (is.na(header)) stop("Cannot locate CellMiner header: ", id)
    headings <- sub("\\s+[a-z]$", "", unlist(raw[header], use.names=FALSE))
    headings[is.na(headings)] <- paste0("unused_", which(is.na(headings)))
    x <- raw[(header+1):.N]; setnames(x, make.unique(headings))
    require_columns(x, "Gene name")
    sample_columns <- grep("^[A-Z]+:", names(x), value=TRUE)
    keys <- normal_key(sub("^[^:]+:", "", sample_columns))
    keys <- sub("ATCC$", "", sub("TB$", "", keys))
    # These source spellings correspond to the explicit NCI60 reference map.
    keys[keys == "MDAMB435"] <- "MDAMB435"
    ids <- sample_lookup$sampleid[match(keys, sample_lookup$sample_key)]
    write_tsv(data.table(source_label=sample_columns, sampleid=ids), paste0(id, "_sample_mapping"))
    keep <- !is.na(ids)
    if (!any(keep)) stop("No CellMiner samples matched")
    matrix <- as.matrix(x[, lapply(.SD, function(v) suppressWarnings(as.numeric(v))), .SDcols=sample_columns[keep]])
    colnames(matrix) <- ids[keep]
    features <- if ("Identifier" %in% names(x)) x$Identifier else paste(x[["Gene name"]], x[["Entrez gene id"]], sep="|")
    profiles[[assay_types[[id]]]] <- gene_se(matrix, features, x[["Gene name"]], assay_types[[id]], units[[id]])
  }
  # DTP includes heterogeneous technologies. Preserve source feature IDs and
  # separate each METHOD / unit combination instead of averaging unlike assays.
  x <- fread(source_file("molecular_target", "\\.TXT$"), na.strings=c("", "NA", "N/A"))
  x[, sampleid := paste(PANELNBR, CELLNBR, sep="|")]
  x <- x[sampleid %in% samples$sampleid]
  x[, value := suppressWarnings(as.numeric(VALUE))]
  x <- x[is.finite(value)]
  groups <- unique(x[, .(ENTITY_MEASURED, METHOD, UNITS)])
  write_tsv(groups, "dtp_assay_groups")
  for (i in seq_len(nrow(groups))) {
    part <- x[groups[i], on=.(ENTITY_MEASURED, METHOD, UNITS)]
    name <- paste0("dtp.", gsub("[^A-Za-z0-9]+", "_", paste(unlist(groups[i]), collapse=".")))
    part[, feature_id := as.character(MOLTID)]
    wide <- dcast(part, feature_id ~ sampleid, value.var="value", fun.aggregate=mean_finite, fill=NA_real_)
    matrix <- as.matrix(wide[, -1]); rownames(matrix) <- wide$feature_id
    info <- unique(part[, .(feature_id, source_gene=GENE, source_title=TITLE)], by="feature_id")
    info <- as.data.frame(info[match(wide$feature_id, feature_id)])
    rownames(info) <- info$feature_id
    profiles[[name]] <- SummarizedExperiment(assays=list(value=matrix), rowData=info,
      colData=S4Vectors::DataFrame(sampleid=colnames(matrix), row.names=colnames(matrix)),
      metadata=list(value_type=groups$UNITS[[i]], feature_type="source molecular-target feature; heterogeneous non-gene targets included"))
  }
}

if (grepl("^GDSC", cfg$dataset$id)) {
  models <- fread(source_file("cmp_models"))
  lookup <- models[!is.na(COSMIC_ID), .(n=uniqueN(model_id), model_id=model_id[[1]]), by=COSMIC_ID][n == 1]
  if (!"model_id" %in% names(samples)) samples[, model_id := NA_character_]
  missing <- is.na(samples$model_id) | samples$model_id == ""
  samples[missing, model_id := lookup$model_id[match(COSMIC_ID, lookup$COSMIC_ID)]]
  lookup <- samples[!is.na(model_id) & model_id != "", .(sampleid, model_id)]
  if (anyDuplicated(lookup$model_id)) stop("Multiple screen samples map to one CMP model")
  cmp_genes <- fread(source_file("cmp_genes"))
  for (path in source_files("rnaseq", "rnaseq_(tpm|fpkm|read_count)_.*\\.csv$")) {
    columns <- names(fread(path, nrows=0, header=TRUE))
    model_columns <- intersect(columns, lookup$model_id)
    x <- fread(path, header=TRUE, select=c("model_id", model_columns))
    x <- x[model_id %in% cmp_genes$gene_id]
    matrix <- as.matrix(x[, lapply(.SD, as.numeric), .SDcols=model_columns])
    colnames(matrix) <- lookup$sampleid[match(model_columns, lookup$model_id)]
    name <- paste0("rnaseq.", sub("rnaseq_(.*)_[0-9]+\\.csv", "\\1", basename(path)))
    profiles[[name]] <- gene_se(matrix, x$model_id, cmp_genes$ensembl_gene_id[match(x$model_id, cmp_genes$gene_id)], name, sub("rnaseq.", "", name, fixed=TRUE))
  }
  x <- fread(source_file("cnv", "WES_pureCN_CNV_genes_[0-9]+\\.csv$"),
    select=c("source", "model_id", "gene_id", "total_copy_number", "seg_mean", "gene_mean", "num_snps", "gatk_mean_log2_copy_ratio"))
  x <- x[source == "Sanger" & model_id %in% lookup$model_id]
  x[, sampleid := lookup$sampleid[match(model_id, lookup$model_id)]]
  for (metric in intersect(c("total_copy_number", "seg_mean", "gene_mean", "num_snps", "gatk_mean_log2_copy_ratio"), names(x))) {
    wide <- dcast(x, gene_id ~ sampleid, value.var=metric, fun.aggregate=mean_finite, fill=NA_real_)
    matrix <- as.matrix(wide[, -1]); name <- paste0("cnv.", metric)
    profiles[[name]] <- gene_se(matrix, wide$gene_id, cmp_genes$ensembl_gene_id[match(wide$gene_id, cmp_genes$gene_id)], name, metric)
  }
  x <- fread(source_file("mutations", "\\.csv$"))
  x <- x[source == "Sanger" & model_id %in% lookup$model_id]
  x[, sampleid := lookup$sampleid[match(model_id, lookup$model_id)]]
  # Retain the variant table; absence of a row is not proof of wild-type.
  write_tsv(x, "gdsc_variants", "data/procdata")
  wide <- dcast(x, gene_symbol ~ sampleid, value.var="gene_symbol", fun.aggregate=length, fill=0L)
  matrix <- as.matrix(wide[, -1]); matrix[matrix>0] <- 1L
  profiles[["mutation.wes"]] <- gene_se(matrix, wide$gene_symbol, wide$gene_symbol, "mutation.wes", "binary reported variant presence; 0 means no reported variant")
  x <- fread(source_file("fusions", "\\.csv$"))
  x <- x[model_id %in% lookup$model_id]
  x[, `:=`(sampleid=lookup$sampleid[match(model_id, lookup$model_id)],
            gene3=gene_ids(cmp_genes$ensembl_gene_id[match(gene_id_3prime, cmp_genes$gene_id)]),
            gene5=gene_ids(cmp_genes$ensembl_gene_id[match(gene_id_5prime, cmp_genes$gene_id)]))]
  write_tsv(x[is.na(gene3) | is.na(gene5)], "fusion_excluded_features")
  x <- x[!is.na(gene3) & !is.na(gene5)]
  x[, feature_id := paste(gene5, gene3, sep="::")]
  wide <- dcast(x, feature_id ~ sampleid, value.var="feature_id", fun.aggregate=function(v) as.integer(length(v)>0), fill=0L)
  matrix <- as.matrix(wide[, -1]); rownames(matrix) <- wide$feature_id
  rd <- as.data.frame(unique(x[, .(feature_id, gene5, gene3)], by="feature_id"))[match(wide$feature_id, unique(x$feature_id)), , drop=FALSE]
  rownames(rd) <- wide$feature_id
  profiles[["fusion.rnaseq"]] <- SummarizedExperiment(assays=list(present=matrix), rowData=rd,
    colData=S4Vectors::DataFrame(sampleid=colnames(matrix), row.names=colnames(matrix)),
    metadata=list(value_type="reported fusion presence", feature_type="5prime::3prime current Ensembl pair"))
}

summary <- rbindlist(lapply(names(profiles), function(n) data.table(assay=n, features=nrow(profiles[[n]]), samples=ncol(profiles[[n]]))))
write_tsv(summary, "molecular_summary")
save_rds(profiles, "data/procdata/molecular_profiles.rds")
