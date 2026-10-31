# Observed and fitted references remain separate. Positive excess denotes
# lower combination viability than the stated reference model.
source("workflow/scripts/common.R")
suppressPackageStartupMessages(library(PharmacoGx))
bundle <- readRDS("data/procdata/curve_response.rds")
x <- bundle$response
mono <- bundle$monotherapy_profiles
keys <- bundle$mono_keys
mono[, four_parameter_dose_support := n_doses >= 4L]
bundle$monotherapy_profiles <- mono
write_tsv(mono, "monotherapy_fits")
# Nikta retains observed NOTZ mono means above one; only curve-fitting inputs are truncated.
parameters <- c("Hill_slope", "E0", "E_inf", "EC50")
for (side in 1:2) {
  part <- copy(mono[, c(keys, parameters), with=FALSE])
  if (side == 2L) setnames(part, "treatment1id", "treatment2id")
  join <- if (side == 1L) keys else sub("treatment1id", "treatment2id", keys, fixed=TRUE)
  renamed <- paste0(parameters, "_", side)
  setnames(part, parameters, renamed)
  x[part, on=join, (renamed) := mget(paste0("i.", renamed))]
  predicted <- paste0("viability_", side, "_fitted")
  dose <- x[[paste0("treatment", side, "dose")]]
  pred <- hillCurve(log10(dose), x[[renamed[1]]], log10(x[[renamed[4]]]),
                    x[[renamed[3]]], E_ninf=x[[renamed[2]]])
  pred[!is.finite(pred)] <- NA_real_
  x[, (predicted) := pmin(1, pmax(0, pred))]
}
x[, `:=`(HSA_observed_expected=pmin(viability_1, viability_2),
          Bliss_observed_expected=viability_1 * viability_2,
          HSA_fitted_expected=pmin(viability_1_fitted, viability_2_fitted),
          Bliss_fitted_expected=viability_1_fitted * viability_2_fitted)]
# These two package APIs assume E0=1, unlike PR178's fitted E0. Their names
# explicitly retain that assumption for comparison with the peer scripts.
x[, `:=`(Loewe_E0_1_expected=NA_real_, ZIP_E0_1_expected=NA_real_)]
if (cfg$dataset$id %in% c("GDSC_MATRIX", "NCI_ALMANAC")) {
  eligible <- which(!is.na(x$treatment2id) & is.finite(x$Hill_slope_1) & is.finite(x$Hill_slope_2) &
    is.finite(x$EC50_1) & is.finite(x$EC50_2) &
    x$EC50_1 > 0 & x$EC50_2 > 0 & x$treatment1dose > 0 & x$treatment2dose > 0 &
    is.finite(x$E_inf_1) & is.finite(x$E_inf_2))
  for (start in seq.int(1L, max(1L, length(eligible)), by=10000L)) {
    if (!length(eligible)) break
    idx <- eligible[start:min(start+9999L, length(eligible))]
    part <- x[idx]
    safe_expected <- function(fun, z) {
      args <- list(treatment1dose=z$treatment1dose, HS_1=z$Hill_slope_1,
        E_inf_1=z$E_inf_1, EC50_1=z$EC50_1, treatment2dose=z$treatment2dose,
        HS_2=z$Hill_slope_2, E_inf_2=z$E_inf_2, EC50_2=z$EC50_2)
      result <- tryCatch(do.call(fun, args), error=function(e) NULL)
      if (is.null(result)) {
        result <- vapply(seq_len(nrow(z)), function(i) tryCatch(
          as.numeric(do.call(fun, lapply(args, `[`, i)))[1],
          error=function(e) NA_real_), numeric(1))
      }
      result[!is.finite(result)] <- NA_real_
      result
    }
    loewe <- safe_expected(computeLoewe, part)
    zip <- safe_expected(computeZIP, part)
    x[idx, `:=`(Loewe_E0_1_expected=loewe, ZIP_E0_1_expected=zip)]
    message("Synergy: ", min(start+9999L, length(eligible)), "/", length(eligible), " eligible wells")
  }
}
expected <- c("HSA_observed", "Bliss_observed", "HSA_fitted", "Bliss_fitted", "Loewe_E0_1", "ZIP_E0_1")
for (name in expected) {
  x[, (paste0(name, "_excess")) := get(paste0(name, "_expected")) - viability]
}
combo_keys <- c("sampleid", "tech_rep", "plateid", "treatment1id", "treatment2id")
if (cfg$dataset$id == "NCI_ALMANAC") combo_keys <- setdiff(combo_keys, "plateid")
if (cfg$dataset$id == "GDSC_MATRIX") combo_keys <- c(combo_keys, "matrixid")
if (cfg$dataset$id == "GDSC_ANCHOR") combo_keys <- c(combo_keys, "treatment2dose", "anchor_doses")
scores <- paste0(expected, "_excess")
combo <- x[!is.na(treatment2id), c(list(mean_viability=mean_finite(viability), n_wells=.N,
    n_dose1=uniqueN(treatment1dose), n_dose2=uniqueN(treatment2dose),
    n_dose_pairs=uniqueN(paste(treatment1dose, treatment2dose)),
    n_missing_observed=sum(!is.finite(viability_1) | !is.finite(viability_2)),
    n_missing_fitted=sum(!is.finite(viability_1_fitted) | !is.finite(viability_2_fitted))),
    setNames(lapply(.SD, mean_finite), paste0("mean_", scores))), by=combo_keys, .SDcols=scores]
if (cfg$dataset$id == "GDSC_MATRIX") combo[, complete_7x7 := n_dose_pairs == 49 & n_dose1 == 7 & n_dose2 == 7]
if (nrow(bundle$combination_fits)) combo <- merge(combo, bundle$combination_fits, by=combo_keys, all.x=TRUE)
write_tsv(combo, "combination_profiles")
bundle$response <- x
bundle$combination_profiles <- combo
save_rds(bundle, "data/procdata/fitted_response.rds")
