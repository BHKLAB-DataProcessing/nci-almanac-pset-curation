# Naming-only adapter for the PR178 PharmacoSet2 / indexed TRE schema.
# Intentionally avoid legacy updateSampleId/updateTreatmentId, which assume
# one sample per column / a singular treatmentid. No constructor or fit runs.
map_ids <- function(x, mapping) {
  y <- unname(mapping[as.character(x)])
  if (any(!is.na(x) & is.na(y))) stop("Unmapped identifier: ", paste(unique(x[!is.na(x) & is.na(y)]),collapse=", "))
  y
}
validate_name_maps <- function(drugs, samples) {
  for (x in list(drugs$nsc, drugs$named_treatmentid, samples$source_sampleid, samples$named_sampleid))
    if (anyNA(x) || any(!nzchar(x)) || anyDuplicated(x)) stop("Naming map must be a complete bijection")
  invisible(TRUE)
}
name_almanac_pset <- function(p, drugs, samples, provenance=list()) {
  validate_name_maps(drugs,samples)
  dm <- setNames(drugs$named_treatmentid, drugs$nsc)
  sm <- setNames(samples$named_sampleid, samples$source_sampleid)
  stopifnot(setequal(p@treatment$treatmentid,names(dm)),setequal(p@sample$sampleid,names(sm)))
  t <- p@treatmentResponse
  stopifnot(is(t,"TreatmentResponseExperiment"), all(c("treatment1id","treatment2id","NSC1","NSC2") %in% names(t@rowData)))
  # These small/index metadata copies protect the input object's data.tables.
  rd <- data.table::copy(t@rowData); cd <- data.table::copy(t@colData)
  for (n in c("treatment1id","treatment2id")) data.table::set(rd,j=n,value=map_ids(rd[[n]],dm))
  data.table::set(cd,j="sampleid",value=map_ids(cd$sampleid,sm))
  rkeys <- as.character(t@.intern$rowIDs); ckeys <- as.character(t@.intern$colIDs)
  rd[, .rownames := do.call(paste,c(.SD,sep=":")), .SDcols=rkeys]
  cd[, .colnames := do.call(paste,c(.SD,sep=":")), .SDcols=ckeys]
  stopifnot(!anyDuplicated(rd,by=rkeys),!anyDuplicated(cd,by=ckeys),!anyDuplicated(rd$.rownames),!anyDuplicated(cd$.colnames))
  t@rowData <- rd; t@colData <- cd
  # In this schema all scientific measurements reside in the four assay
  # tables; IDs reside exclusively in rowData/colData. Fail on schema drift.
  stopifnot(!any(vapply(t@assays,function(a)any(c("sampleid","treatment1id","treatment2id") %in% names(a)),logical(1))))
  p@treatmentResponse <- t
  sd <- p@sample; dd <- p@treatment
  sd$sampleid_source <- sd$sampleid; sd$sample_name_source <- sd$sample_name
  si <- match(sd$sampleid,samples$source_sampleid)
  sd$NCI60_2026.sampleid <- samples$nci60_2026_sampleid[si]
  sd$sampleid <- map_ids(sd$sampleid,sm); sd$sample_name <- sd$sampleid
  rownames(sd) <- sd$sampleid
  dd$NSC_number <- dd$treatmentid; dd$treatmentid_source <- dd$treatmentid
  dd$treatment_name_source <- dd$treatment_name
  di <- match(dd$treatmentid,drugs$nsc)
  dd$treatmentid <- map_ids(dd$treatmentid,dm); dd$treatment_name <- dd$treatmentid
  dd$source_formulation_name <- drugs$preferred_name[di]
  dd$source_CAS <- drugs$release_cas[di]
  dd$identity_review_flags <- drugs$identity_review_flags[di]
  rownames(dd) <- dd$treatmentid
  p@sample <- sd; p@treatment <- dd
  for (kind in c("sample","treatment")) {
    field <- paste0(kind,"id"); m <- if(kind=="sample") sm else dm
    z <- p@curation[[kind]]
    z[[paste0(field,"_source")]] <- z[[field]]
    z[[field]] <- map_ids(z[[field]],m); rownames(z) <- z[[field]]
    p@curation[[kind]] <- z
  }
  mae <- p@molecularProfiles
  exps <- as.list(MultiAssayExperiment::experiments(mae))
  for (n in names(exps)) {
    se <- exps[[n]]
    old <- colnames(se); colnames(se) <- map_ids(old,sm)
    SummarizedExperiment::colData(se)$sampleid_source <- old
    if ("sampleid" %in% names(SummarizedExperiment::colData(se)))
      SummarizedExperiment::colData(se)$sampleid <- map_ids(SummarizedExperiment::colData(se)$sampleid,sm)
    exps[[n]] <- se
  }
  mc <- mae@colData
  mc$sampleid_source <- as.character(mc$sampleid)
  mc$sample_name_source <- as.character(mc$sample_name)
  mc$sampleid <- map_ids(mc$sampleid,sm); mc$sample_name <- mc$sampleid
  rownames(mc) <- map_ids(rownames(mc),sm)
  mm <- mae@sampleMap
  mm$primary <- map_ids(mm$primary,sm); mm$colname <- map_ids(mm$colname,sm)
  mae@ExperimentList <- MultiAssayExperiment::ExperimentList(exps)
  mae@colData <- mc; mae@sampleMap <- mm
  p@molecularProfiles <- mae
  p@annotation$naming_revision <- c(provenance,list(
    scope="Naming only; all scientific values, indexed assay membership and source NSCs unchanged; no refit or potency patch",
    treatment_map=as.data.frame(drugs),sample_map=as.data.frame(samples)))
  methods::validObject(p);methods::validObject(t);methods::validObject(mae)
  p
}
scientific_fingerprint <- function(p) {
  hash <- function(v) digest::digest(v,algo="sha256",serialize=TRUE)
  table_hash <- function(d,omit=character()) lapply(as.list(d)[setdiff(names(d),omit)],hash)
  t <- p@treatmentResponse
  exps <- as.list(MultiAssayExperiment::experiments(p@molecularProfiles))
  list(row=table_hash(t@rowData,c("treatment1id","treatment2id",".rownames")),
       col=table_hash(t@colData,c("sampleid",".colnames")),
       assays=lapply(t@assays,table_hash), index=table_hash(t@.intern$assayIndex),
       assay_keys=lapply(t@.intern$assayKeys,as.character),
       molecular=lapply(exps,function(se)list(
         values=lapply(as.list(SummarizedExperiment::assays(se,withDimnames=FALSE)),hash),
         row_data=hash(SummarizedExperiment::rowData(se)),row_names=hash(rownames(se)))) )
}
# Select both orientations of an exact named pair through TRE's supported
# row/column index interface. No singleton treatment assumptions.
named_pair_response <- function(p,drug1,drug2,samples=NULL) {
  stopifnot(all(c(drug1,drug2) %in% PharmacoGx::drugNames(p)))
  t <- p@treatmentResponse; rd <- t@rowData; cd <- t@colData
  ri <- which((rd$treatment1id==drug1 & rd$treatment2id==drug2) |
              (rd$treatment1id==drug2 & rd$treatment2id==drug1))
  ck <- unique(t@.intern$assayIndex$colKey[t@.intern$assayIndex$rowKey %in% rd$rowKey[ri]])
  ci <- which(cd$colKey %in% ck)
  if(!is.null(samples)) {stopifnot(all(samples %in% PharmacoGx::sampleNames(p)));ci <- ci[cd$sampleid[ci] %in% samples]}
  t[ri,ci]
}
