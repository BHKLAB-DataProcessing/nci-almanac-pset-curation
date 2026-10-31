.libPaths(.Library)
suppressPackageStartupMessages(library(data.table))
cfg <- yaml::read_yaml("config/pipeline.yaml")
source("workflow/scripts/check_environment.R")
check_curation_environment(cfg$packages)
data.table::setDTthreads(1L)
members <- rbindlist(lapply(list.files("data/procdata/extracted", "^manifest.tsv$", recursive=TRUE, full.names=TRUE), fread))
source_files <- function(id, pattern = NULL) {
  paths <- members[source_id == id, path]
  if (!is.null(pattern)) paths <- paths[grepl(pattern, paths)]
  if (!length(paths)) stop("No source member for ", id, " / ", pattern)
  paths
}
source_file <- function(id, pattern = NULL) {
  paths <- source_files(id, pattern)
  if (length(paths) != 1L) stop("Ambiguous source member for ", id, ": ", paste(paths, collapse=", "))
  paths[[1]]
}
require_columns <- function(x, columns) {
  absent <- setdiff(columns, names(x))
  if (length(absent)) stop("Missing required source columns: ", paste(absent, collapse=", "))
}
save_rds <- function(x, path) {
  dir.create(dirname(path), recursive=TRUE, showWarnings=FALSE)
  temporary <- paste0(path, ".part")
  saveRDS(x, temporary, compress="gzip")
  if (!file.rename(temporary, path)) stop("Cannot replace ", path)
}
write_tsv <- function(x, name, directory="data/procdata/metadata") {
  dir.create(directory, recursive=TRUE, showWarnings=FALSE)
  fwrite(x, file.path(directory, paste0(name, ".tsv")), sep="\t", na="NA", quote="auto")
}
normal_key <- function(x) toupper(gsub("[^A-Za-z0-9]", "", x))
mean_finite <- function(x) if (any(is.finite(x))) mean(x[is.finite(x)]) else NA_real_
collapse_text <- function(x) paste(sort(unique(x[!is.na(x) & nzchar(x)])), collapse="|")
canonical_response <- function(x) {
  require_columns(x, c("sampleid", "treatment1id", "treatment1dose", "viability", "plateid"))
  if (!"treatment2id" %in% names(x)) x[, treatment2id := NA_character_]
  if (!"treatment2dose" %in% names(x)) x[, treatment2dose := NA_real_]
  if (!"tech_rep" %in% names(x)) x[, tech_rep := plateid]
  if (!"wellid" %in% names(x)) x[, wellid := as.character(seq_len(.N)), by=.(sampleid, plateid)]
  for (column in c("sampleid", "treatment1id", "treatment2id", "plateid", "tech_rep", "wellid"))
    set(x, j=column, value=as.character(x[[column]]))
  x[, responseid := as.character(seq_len(.N))]
  if (anyNA(x$sampleid) || anyNA(x$treatment1id)) stop("Missing sample or primary treatment IDs")
  x
}
