# Current gene-level identifiers are constrained to the configured GENCODE release.
gtf <- fread(source_file("gencode"), sep="\t", header=FALSE, quote="",
             select=c(1L,3L,4L,5L,7L,9L))
setnames(gtf, c("chromosome", "feature", "start", "end", "strand", "attributes"))
genes <- gtf[feature == "gene", .(gene_id=sub('.*gene_id "([^"]+)".*', '\\1', attributes),
                                  gene_primary_name=sub('.*gene_name "([^"]+)".*', '\\1', attributes),
                                  chromosome, start, end, strand)]
genes[, gene_id := sub("\\.[0-9]+$", "", gene_id)]
genes <- unique(genes, by="gene_id")
rm(gtf); gc(FALSE)
symbol_lookup <- genes[, .(n_ids=uniqueN(gene_id), gene_id=gene_id[[1]]), by=gene_primary_name][n_ids == 1]
gene_ids <- function(ids) {
  ids <- sub("\\.[0-9]+$", "", as.character(ids))
  mapped <- ids
  idx <- !grepl("^ENSG[0-9]+$", ids) & !is.na(ids)
  mapped[idx] <- symbol_lookup$gene_id[match(ids[idx], symbol_lookup$gene_primary_name)]
  mapped[!mapped %in% genes$gene_id] <- NA_character_
  mapped
}
gene_se <- function(matrix, source_features, gene_candidates, name, value_type) {
  stopifnot(nrow(matrix) == length(source_features), length(source_features) == length(gene_candidates))
  mapping <- data.table(source_feature=as.character(source_features), gene_id=gene_ids(gene_candidates),
                        source_index=seq_along(source_features))
  setorder(mapping, gene_id, source_feature, source_index)
  mapping[, exclusion_reason := fifelse(is.na(gene_id), "no_unique_current_Ensembl_match", "")]
  mapping[!is.na(gene_id), exclusion_reason := fifelse(duplicated(gene_id), "duplicate_gene_keep_lexical_first_feature", "")]
  write_tsv(mapping, paste0(name, "_feature_mapping"))
  gene_audit <- copy(mapping)
  gene_audit[, `:=`(In.AnnotationDB=!is.na(gene_id), PubChem.CID=NA_character_,
                    AnnotationGx.Mapped=NA, AnnotationGx.Title=genes$gene_primary_name[match(gene_id, genes$gene_id)])]
  write_tsv(gene_audit[In.AnnotationDB == FALSE], paste0(name, "_annotationdb_gene_misses"))
  keep <- mapping[exclusion_reason == ""]
  if (!nrow(keep) || !ncol(matrix)) stop("No mapped molecular features/samples for ", name)
  matrix <- matrix[keep$source_index, , drop=FALSE]
  rownames(matrix) <- keep$gene_id
  rd <- as.data.frame(merge(keep[, .(gene_id, source_feature)], genes, by="gene_id", sort=FALSE))
  rd <- rd[match(rownames(matrix), rd$gene_id), , drop=FALSE]
  rownames(rd) <- rd$gene_id
  SummarizedExperiment::SummarizedExperiment(assays=list(value=matrix), rowData=rd,
    colData=S4Vectors::DataFrame(sampleid=colnames(matrix), row.names=colnames(matrix)),
    metadata=list(value_type=value_type))
}
