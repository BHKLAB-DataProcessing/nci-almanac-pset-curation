source("workflow/scripts/common.R")
bundle <- readRDS("data/procdata/response.rds")
samples <- copy(bundle$sample)
drugs <- copy(bundle$treatment)
source("workflow/scripts/annotation_cache.R")
cell_cache <- source_file("cellosaurus")
cell_url <- paste0(sub("/$","",cfg$annotationdb$base_url),"/cell_line/all")
annotation_receipts <- list(annotationdb_read_response(cell_url,cell_cache,"cell_line_all",cfg))
cells <- as.data.table(jsonlite::fromJSON(cell_cache))
cells[, key := normal_key(name)]
cells <- cells[, .(n_ids=uniqueN(accession), accession=accession[[1]], title=name[[1]]), by=key][n_ids == 1]
sample_keys <- samples$sample_name
if (cfg$dataset$id == "HMCL")
  sample_keys <- sub("_(DSMZ|ATCC|PLB|JCRB|ECACC|RIKEN|DJ|JJKscc[A-Z0-9]+|sus)$", "", sample_keys)
idx <- match(normal_key(sample_keys), cells$key)
samples[, `:=`(Cellosaurus.Accession=cells$accession[idx], In.AnnotationDB=!is.na(idx),
                 AnnotationGx.Mapped=NA, AnnotationGx.Title=cells$title[idx])]
# Public CMP identifiers connect the combination screen to independent omics.
if (grepl("^GDSC", cfg$dataset$id)) {
  models <- fread(source_file("cmp_models"))
  models[, cosmic_key := as.character(COSMIC_ID)]
  lookup <- models[!is.na(cosmic_key) & cosmic_key != "", .(n=uniqueN(model_id), model_id=model_id[[1]], tissue=tissue[[1]]), by=cosmic_key][n == 1]
  if (!"model_id" %in% names(samples)) samples[, model_id := NA_character_]
  j <- match(samples$COSMIC_ID, lookup$cosmic_key)
  missing <- is.na(samples$model_id) | samples$model_id == ""
  samples[missing, model_id := lookup$model_id[j[missing]]]
  if (!"tissueid" %in% names(samples)) samples[, tissueid := lookup$tissue[j]]
  write_tsv(samples[is.na(model_id) | model_id == ""], "unmapped_cmp_samples")
}
drugs[, `:=`(PubChem.CID=NA_character_, In.AnnotationDB=FALSE, AnnotationGx.Mapped=NA,
               AnnotationGx.Title=NA_character_)]
dir.create("data/rawdata/metadata/annotationdb", recursive=TRUE, showWarnings=FALSE)
queries <- sort(unique(drugs$treatment_name[!is.na(drugs$treatment_name) &
                         !grepl("[|+]", drugs$treatment_name)]))
candidates <- list()
# Batch requests use the API's documented repeated `compound` parameter.
# Only exact returned titles / mapped identifiers are accepted, never row order.
for (start in seq.int(1L, max(1L, length(queries)), by=20L)) {
  batch <- queries[start:min(start+19L, length(queries))]
  batch <- batch[!is.na(batch)]
  if (!length(batch)) next
  key <- digest::digest(batch, algo="sha256")
  cache <- file.path("data/rawdata/metadata/annotationdb", paste0(key, ".json"))
  url <- paste0(sub("/$", "", cfg$annotationdb$base_url), "/compound/many?",
                paste(paste0("compound=", vapply(batch, utils::URLencode, "", reserved=TRUE)), collapse="&"))
  annotation_receipts[[length(annotation_receipts)+1L]] <- annotationdb_read_response(url,cache,key,cfg)
  records <- jsonlite::fromJSON(cache, simplifyVector=FALSE)
  for (record in records) {
    identifiers <- c(record$title, strsplit(record$mapped_name, ";\\s*")[[1]])
    hit <- which(normal_key(drugs$treatment_name) %in% normal_key(identifiers))
    if (length(hit)) {
      candidates[[length(candidates)+1L]] <- data.table(treatmentid=drugs$treatmentid[hit],
                            PubChem.CID=as.character(record$cid), AnnotationGx.Title=record$title)
    }
  }
}
if (length(candidates)) {
  candidate_table <- unique(rbindlist(candidates))
  write_tsv(candidate_table, "compound_mapping_candidates")
  unambiguous <- candidate_table[, .(n_candidates=uniqueN(PubChem.CID),
      PubChem.CID=PubChem.CID[[1]], AnnotationGx.Title=AnnotationGx.Title[[1]]), by=treatmentid][n_candidates == 1]
  drugs[unambiguous, on="treatmentid", `:=`(PubChem.CID=i.PubChem.CID, AnnotationGx.Title=i.AnnotationGx.Title, In.AnnotationDB=TRUE)]
}
write_tsv(samples, "sample_annotations")
write_tsv(drugs, "treatment_annotations")
write_tsv(samples[In.AnnotationDB == FALSE], "annotationdb_cell_line_misses")
write_tsv(drugs[In.AnnotationDB == FALSE], "annotationdb_compound_misses")
write_tsv(rbindlist(annotation_receipts), "annotationdb_responses")
save_rds(list(sample=samples, treatment=drugs), "data/procdata/annotations.rds")
