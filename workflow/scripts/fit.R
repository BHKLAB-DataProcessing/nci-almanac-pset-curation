# Fits use the same PharmacoGx Hill/IC50/AAC path as the peer scripts.
source("workflow/scripts/common.R")
suppressPackageStartupMessages(library(PharmacoGx))
bundle <- readRDS("data/procdata/response.rds")
x <- bundle$response
min_points <- cfg$parameters$fitting$min_points
fit_one <- function(dose, response) {
  keep <- is.finite(dose) & dose > 0 & is.finite(response)
  dose <- dose[keep]; response <- response[keep]
  empty <- list(Hill_slope=NA_real_, E0=NA_real_, Rsquare=NA_real_, E_inf=NA_real_, EC50=NA_real_, IC50=NA_real_, AUC=NA_real_, AAC=NA_real_,
                n_doses=uniqueN(dose), fit_status="insufficient_doses")
  if (uniqueN(dose) < min_points) return(empty)
  tryCatch({
    fit <- logLogisticRegression(dose, response, viability_as_pct=FALSE, verbose=FALSE)
    list(Hill_slope=as.numeric(fit[["HS"]]), E0=as.numeric(fit[["E0"]]),
      Rsquare=if (length(attr(fit, "Rsquare"))) as.numeric(attr(fit, "Rsquare")) else NA_real_, E_inf=as.numeric(fit[["E_inf"]]), EC50=as.numeric(fit[["EC50"]]),
      IC50=as.numeric(computeIC50(dose, Hill_fit=fit, viability_as_pct=FALSE, verbose=FALSE)),
      AUC=as.numeric(suppressMessages(computeAUC(dose, Hill_fit=fit, viability_as_pct=FALSE, verbose=FALSE))),
      AAC=as.numeric(computeAAC(dose, Hill_fit=fit, viability_as_pct=FALSE, verbose=FALSE)),
      n_doses=uniqueN(dose), fit_status="ok")
  }, error=function(e) {empty$fit_status <- paste0("error: ", conditionMessage(e)); empty})
}
fit_groups <- function(data, keys, stage) {
  if (!nrow(data)) return(data.table())
  signature <- digest::digest(list(input=digest::digest(file="data/procdata/response.rds", algo="sha256"),
    settings=cfg$parameters$fitting, package=utils::packageDescription("PharmacoGx")$RemoteSha,
    fitter=digest::digest(fit_one, algo="sha256"), keys=keys,
    data=digest::digest(data, algo="sha256"), stage=stage), algo="sha256")
  directory <- file.path("data/procdata/fit_checkpoints", signature)
  dir.create(directory, recursive=TRUE, showWarnings=FALSE)
  data <- copy(data)
  data[, fit_chunk := (.GRP-1L) %/% cfg$parameters$fitting$chunk_size, by=keys]
  chunks <- sort(unique(data$fit_chunk))
  outputs <- parallel::mclapply(chunks, function(chunk) {
    path <- file.path(directory, paste0(chunk, ".rds"))
    if (!file.exists(path)) {
      result <- data[fit_chunk == chunk, fit_one(treatment1dose, pmin(1, pmax(0, viability))), by=keys]
      save_rds(result, path)
    }
    message(stage, ": completed chunk ", chunk+1L, "/", length(chunks))
    path
  }, mc.cores=as.integer(cfg$parameters$fitting$workers), mc.preschedule=FALSE)
  if (any(vapply(outputs, inherits, logical(1), "try-error"))) stop("Curve-fitting worker failed")
  rbindlist(lapply(outputs, readRDS), use.names=TRUE, fill=TRUE)
}
# Average technical wells within the same source context before fitting.
mono <- x[is.na(treatment2id), .(viability=mean_finite(viability), n_wells=.N),
          by=.(sampleid, tech_rep, plateid, treatment1id, treatment1dose)]
if (cfg$dataset$id == "GDSC_MATRIX") {
  # Source mono references repeat across combination wells. Deduplicate by
  # plate/drug/dose/observed value so combination coverage cannot reweight them.
  mono <- unique(rbind(x[, .(sampleid, tech_rep, plateid, matrixid, treatment1id, treatment1dose, viability=viability_1)],
                      x[, .(sampleid, tech_rep, plateid, matrixid, treatment1id=treatment2id,
                            treatment1dose=treatment2dose, viability=viability_2)]))
  mono <- mono[, .(viability=mean_finite(viability), n_wells=.N),
               by=.(sampleid, tech_rep, plateid, matrixid, treatment1id, treatment1dose)]
}
mono_keys <- c("sampleid", "tech_rep", "plateid", "treatment1id")
if (cfg$dataset$id == "GDSC_MATRIX") mono_keys <- c(mono_keys, "matrixid")
if (cfg$dataset$id == "NCI_ALMANAC") {
  # The reference pools monotherapy within screener/study, because the
  # matching single-agent curves can be located on different plates.
  mono <- x[is.na(treatment2id), .(viability=mean_finite(viability), n_wells=.N),
            by=.(sampleid, tech_rep, treatment1id, treatment1dose)]
  mono_keys <- setdiff(mono_keys, "plateid")
}
mono_profiles <- fit_groups(mono, mono_keys, "monotherapy")
write_tsv(mono_profiles, "monotherapy_fits")
if (!"viability_1" %in% names(x)) {
  x[mono, on=c(mono_keys, "treatment1dose"), viability_1 := i.viability]
}
if (!"viability_2" %in% names(x)) {
  join <- c(setNames(setdiff(mono_keys, "treatment1id"), setdiff(mono_keys, "treatment1id")),
    treatment2id="treatment1id", treatment2dose="treatment1dose")
  x[mono, on=join, viability_2 := i.viability]
}
bundle$response <- x
bundle$monotherapy_profiles <- mono_profiles
bundle$mono_keys <- mono_keys
bundle$combination_fits <- if (cfg$dataset$id == "GDSC_ANCHOR") {
  fit_groups(x[!is.na(treatment2id)],
    c("sampleid", "tech_rep", "plateid", "treatment1id", "treatment2id", "treatment2dose", "anchor_doses"), "combination")
} else data.table()
if (nrow(mono_profiles) && !any(mono_profiles$fit_status == "ok")) stop("No successful monotherapy fits")
save_rds(bundle, "data/procdata/curve_response.rds")
