# Reproduce Nikta nci_almanac_curation_optimised.R (2026-10-23) from public raw measurements.
# The peer PSet is never a production input; original NSCs remain primary IDs.
source("workflow/scripts/common.R")
x <- fread(source_file("screen", "\\.csv$"), colClasses=c(NSC1="character", NSC2="character"))
require_columns(x, c("NSC1", "NSC2", "CONC1", "CONC2", "PERCENTGROWTH", "PERCENTGROWTHNOTZ", "CELLNAME", "SCREENER", "STUDY", "PLATE"))
if (any(!is.na(x$CONCUNIT1) & !x$CONCUNIT1 %in% c("", "M")) ||
    any(!is.na(x$CONCUNIT2) & !x$CONCUNIT2 %in% c("", "M"))) stop("Unexpected ALMANAC concentration unit")
samples <- unique(x[, .(sampleid=paste(PANELNBR, CELLNBR, sep="|"), sample_name=CELLNAME, tissueid=PANEL)])
if (anyDuplicated(samples$sampleid)) stop("Conflicting NCI panel/cell IDs")
long <- x[, .(sampleid=paste(PANELNBR, CELLNBR, sep="|"), treatment1id=NSC1, treatment2id=NSC2,
  treatment1dose=CONC1, treatment2dose=CONC2, dose_unit="M",
  NSC1=NSC1, NSC2=NSC2,
  PERCENTGROWTH=PERCENTGROWTH, PERCENTGROWTHNOTZ=PERCENTGROWTHNOTZ,
  viability_reported=PERCENTGROWTHNOTZ / 100, viability=pmax(0, PERCENTGROWTHNOTZ / 100),
  TESTVALUE=TESTVALUE, CONTROLVALUE=CONTROLVALUE, TZVALUE=TZVALUE,
  expected_growth_published=EXPECTEDGROWTH, score_published=SCORE,
  valid_source=VALID, plateid=as.character(PLATE), tech_rep=paste(SCREENER, STUDY, sep="_"),
  wellid=as.character(COMBODRUGSEQ))]
long[treatment2id == "", `:=`(treatment2id=NA_character_, treatment2dose=NA_real_)]
long[, `:=`(is_mono=is.na(treatment2id) | is.na(treatment2dose))]
long[, is_combo := !is_mono]
long <- canonical_response(long)
source("workflow/scripts/parse_compound_names.R")
drugnames <- parse_compound_names(source_file("combocompoundnames_all"))
ids <- unique(c(long$treatment1id, na.omit(long$treatment2id)))
drugs <- merge(data.table(treatmentid=ids), unique(drugnames, by="treatmentid"), by="treatmentid", all.x=TRUE)
# First source synonym remains the baseline display field. The explicit,
# reviewed naming policy is applied during final curation.
if (anyNA(drugs$treatment_name)) stop("Missing public ALMANAC compound names")
write_tsv(drugnames, "source_treatment_aliases")
write_tsv(x[, .N, by=VALID], "almanac_source_validity")
save_rds(list(response=long, sample=samples, treatment=drugs), "data/procdata/response.rds")
