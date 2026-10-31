# Apply only the reviewed canonical-treatment and RNA-seq interface changes.
# All indexed response tables, source NSC keys and molecular matrices stay intact.
almanac_compatibility_v2 <- function(p, provenance=list()) {
  old_treatments <- p@treatment
  old_curation <- p@curation$treatment
  vem <- c("Vemurafenib [NSC 753082]","Vemurafenib [NSC 761431]")
  stopifnot(all(vem %in% old_treatments$treatmentid),!"Vemurafenib" %in% old_treatments$treatmentid)
  rename <- function(x) {x[x %in% vem] <- "Vemurafenib";x}
  t <- p@treatmentResponse; rd <- data.table::copy(t@rowData)
  for(n in c("treatment1id","treatment2id")) data.table::set(rd,j=n,value=rename(rd[[n]]))
  # NSC1/NSC2 remain in rowIDs, so canonical naming does not collapse source rows.
  stopifnot(all(c("NSC1","NSC2") %in% as.character(t@.intern$rowIDs)))
  rd[, .rownames := do.call(paste,c(.SD,sep=":")), .SDcols=as.character(t@.intern$rowIDs)]
  stopifnot(!anyDuplicated(rd,by=as.character(t@.intern$rowIDs)),!anyDuplicated(rd$.rownames))
  t@rowData <- rd;p@treatmentResponse <- t
  combine <- function(d,key) {
    groups <- split(seq_len(nrow(d)),factor(d[[key]],levels=unique(d[[key]])))
    ans <- lapply(groups,function(ii) {
      row <- d[ii[1],,drop=FALSE]
      for(n in names(d)) {
        values <- unique(d[[n]][ii])
        if(length(values)>1L) {
          if(!is.character(values)) stop("Conflicting non-character annotation: ",n)
          row[[n]] <- paste(sort(values[!is.na(values) & nzchar(values)]),collapse="///")
        }
      }
      row
    })
    ans <- do.call(rbind,ans);rownames(ans)<-ans[[key]];ans
  }
  d <- old_treatments
  d$treatmentid <- rename(d$treatmentid);d$treatment_name<-d$treatmentid
  d$name <- d$treatmentid;d$treatmentid_dataset<-as.character(d$NSC_number)
  p@treatment <- combine(d,"treatmentid")
  c <- old_curation;c$treatmentid <- rename(c$treatmentid)
  p@curation$treatment <- combine(c,"treatmentid")
  stopifnot(p@treatment["Vemurafenib","treatmentid_dataset"]=="753082///761431",
            p@curation$treatment["Vemurafenib","treatmentid_dataset"]=="753082///761431")
  mae <- p@molecularProfiles
  stopifnot("rnaseq.composite" %in% names(mae),!"rnaseq.comp" %in% names(mae))
  names(mae)[names(mae)=="rnaseq.composite"] <- "rnaseq.comp"
  se <- mae[["rnaseq.comp"]]
  stopifnot(identical(SummarizedExperiment::assayNames(se),"value"))
  S4Vectors::metadata(se)$source_experiment_name <- "rnaseq.composite"
  S4Vectors::metadata(se)$source_inner_assay_name <- "value"
  S4Vectors::metadata(se)$annotation <- "rnaseq.comp"
  SummarizedExperiment::assayNames(se) <- "exprs"
  mae[["rnaseq.comp"]] <- se;p@molecularProfiles <- mae
  p@annotation$compatibility_v2 <- c(provenance,list(
    changes=c("Canonical Vemurafenib metadata with source IDs 753082///761431; no response aggregation",
      "rnaseq.composite/value renamed rnaseq.comp/exprs; no type or value change"),
    source_treatment_metadata=old_treatments,source_treatment_curation=old_curation,
    limitations=c("No rna microarray or rnaseq.iso data added", "Default sensitivityProfiles accessor remains incompatible with context-level profiles; use explicit sProfiles helper",
      "Existing fitted-potency caveats unchanged; no refit or package patch")))
  methods::validObject(p);methods::validObject(t);methods::validObject(mae)
  p
}

# Normalize only the explicitly renamed experiment/list labels when comparing
# fingerprints. Numeric values, source indices, dose/well keys and feature IDs
# still require exact SHA256 equality.
compatibility_scientific_fingerprint <- function(p) {
  result <- scientific_fingerprint(p)
  i <- which(names(result$molecular)=="rnaseq.comp")
  if(length(i)) {
    names(result$molecular)[i] <- "rnaseq.composite"
    stopifnot(identical(names(result$molecular[[i]]$values),"exprs"))
    names(result$molecular[[i]]$values) <- "value"
  }
  result
}
