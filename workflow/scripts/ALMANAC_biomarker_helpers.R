# Companion for the ALMANAC compatibility-v2 RDS. No file writes or refitting.
# Example:
# library(PharmacoGx)
# source("ALMANAC_biomarker_helpers.R")
# p <- readRDS("NCI_ALMANAC_PSet.rds")
# result <- almanac_pair_signature(p, c("Paclitaxel", "Nilotinib"),
#   features=c("ENSG00000000003","ENSG00000000005","ENSG00000000419"))
# result$signature
# Uses the actual RNA-seq layer, not an assumed microarray experiment named rna.

almanac_pair_profiles <- function(pset, drugs, source_nscs=NULL) {
  stopifnot(length(drugs)==2L,all(drugs %in% PharmacoGx::drugNames(pset)))
  info <- PharmacoGx::treatmentInfo(pset)
  ambiguous <- grepl("///",info[drugs,"treatmentid_dataset"],fixed=TRUE)
  if(any(ambiguous) && is.null(source_nscs))
    stop("This pair includes multiple source NSCs under a canonical name. Supply source_nscs in the same order as drugs; no silent pooling.")
  if(!is.null(source_nscs)) {
    stopifnot(length(source_nscs)==2L,!anyNA(source_nscs))
    for(i in 1:2) if(!as.character(source_nscs[i]) %in% strsplit(info[drugs[i],"treatmentid_dataset"],"///",fixed=TRUE)[[1]])
      stop("source_nscs does not match the selected canonical treatment")
  }
  t <- PharmacoGx::treatmentResponse(pset);rd<-t@rowData;cd<-t@colData
  forward <- rd$treatment1id==drugs[1] & rd$treatment2id==drugs[2]
  reverse <- rd$treatment1id==drugs[2] & rd$treatment2id==drugs[1]
  if(!is.null(source_nscs)) {
    forward <- forward & rd$NSC1==source_nscs[1] & rd$NSC2==source_nscs[2]
    reverse <- reverse & rd$NSC1==source_nscs[2] & rd$NSC2==source_nscs[1]
  }
  ri <- which(forward | reverse)
  if(!length(ri)) stop("No observed source pair matches the requested drugs/NSCs")
  ci <- which(cd$colKey %in% unique(t@.intern$assayIndex$colKey[t@.intern$assayIndex$rowKey %in% rd$rowKey[ri]]))
  selected <- t[ri,ci]
  selected[["profiles"]]
}

almanac_pair_signature <- function(pset,drugs,features,mDataType="rnaseq.comp",
                                   source_nscs=NULL,measure="mean_Bliss_observed_excess") {
  for(pkg in c("PharmacoGx","SummarizedExperiment","data.table"))
    if(!requireNamespace(pkg,quietly=TRUE)) stop("Required package unavailable: ",pkg)
  if(is.na(parallel::detectCores())) stop("CPU detection unavailable in this execution sandbox; run in a normal R session. The package requires detectCores() even with nthread=1.")
  if(!mDataType %in% names(pset@molecularProfiles)) stop("Requested molecular experiment is absent")
  se <- pset@molecularProfiles[[mDataType]]
  if(!all(features %in% rownames(se))) stop("Requested feature IDs are absent; use this source's rownames")
  d <- data.table::as.data.table(almanac_pair_profiles(pset,drugs,source_nscs))
  if(!measure %in% names(d) || !is.numeric(d[[measure]])) stop("Requested profile measure is absent or nonnumeric")
  # Unit of aggregation is the existing fitted/observed technical context.
  # NSC selection happens before equal-weight averaging across contexts.
  good <- is.finite(d[[measure]])
  z <- data.table::data.table(sampleid=d$sampleid[good],value=d[[measure]][good])
  summary <- z[,list(value=mean(value),technical_contexts=.N),by="sampleid"]
  if(!nrow(summary)) stop("No finite response values for the requested pair")
  shared <- intersect(summary$sampleid,colnames(se))
  if(length(shared)<3L) stop("Fewer than three response/molecular shared samples")
  s <- matrix(summary$value,nrow=1,dimnames=list(paste(drugs,collapse=" + "),summary$sampleid))
  sig <- PharmacoGx::drugSensitivitySig(pset,mDataType=mDataType,
    sProfiles=s,features=features,cells=shared,
    modeling.method="pearson",inference.method="analytic",
    molecular.summary.stat="mean",sensitivity.summary.stat="mean",
    standardize="SD",nthread=1,verbose=FALSE)
  list(signature=sig,response_summary=as.data.frame(summary),shared_samples=shared,
       excluded_samples=setdiff(summary$sampleid,shared),
       arguments=list(drugs=drugs,source_nscs=source_nscs,measure=measure,
         mDataType=mDataType,features=features,aggregation="Equal-weight mean across selected technical-context profiles per sample; no fitted potency recomputation",
         modeling.method="pearson",inference.method="analytic",standardize="SD",nthread=1))
}
