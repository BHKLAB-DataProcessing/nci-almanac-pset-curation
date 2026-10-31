# Build complete maps from original NCI source files plus reviewed label choices.
almanac_source_name_maps <- function(samples, treatments,
    names_file=source_file("combocompoundnames_small"),
    structures_file=source_file("combocompoundset"), policy_file="config/naming-policy.json") {
  policy <- jsonlite::fromJSON(policy_file,simplifyVector=FALSE)
  source("workflow/scripts/parse_compound_names.R",local=TRUE)
  names <- parse_compound_names(names_file)
  preferred <- names[!duplicated(treatmentid)]
  ids <- as.character(treatments$treatmentid)
  stopifnot(!anyDuplicated(ids),all(ids %in% preferred$treatmentid))
  drugs <- data.frame(nsc=ids,preferred_name=preferred$treatment_name[match(ids,preferred$treatmentid)],stringsAsFactors=FALSE)
  replace_values <- function(values,ids,overrides) {
    hit <- match(ids,names(overrides));use <- !is.na(hit)
    values[use] <- unlist(overrides[hit[use]],use.names=FALSE);values
  }
  drugs$preferred_name <- replace_values(drugs$preferred_name,ids,policy$formulation_overrides)
  drugs$named_treatmentid <- replace_values(drugs$preferred_name,ids,policy$drug_display_overrides)
  # Source CAS annotations are parsed directly, not copied from a local PSet.
  blocks <- strsplit(paste(readLines(structures_file,warn=FALSE),collapse="\n"),"$$$$",fixed=TRUE)[[1]]
  field <- function(lines,name) {
    i <- which(trimws(lines)==paste0(">  <",name,">"))
    if(!length(i)) return(NA_character_)
    trimws(lines[i[[1]]+1L])
  }
  sdf <- lapply(blocks,function(block) {
    lines <- strsplit(trimws(block),"\n",fixed=TRUE)[[1]]
    data.frame(nsc=field(lines,"NSC"),cas=field(lines,"CAS"))
  })
  sdf <- do.call(rbind,sdf);sdf<-sdf[!is.na(sdf$nsc),,drop=FALSE]
  if(anyDuplicated(sdf$nsc)) stop("Duplicate NCI structure NSCs")
  drugs$release_cas <- sdf$cas[match(ids,sdf$nsc)]
  drugs$identity_review_flags <- ifelse(ids %in% c("753082","761431"),"distinct_source_NSCs_retained_without_structure_equivalence_claim",NA_character_)
  for (id in c("122758","125066")) drugs$identity_review_flags[ids==id] <-
    if(id=="122758") "tretinoin_not_isotretinoin" else "bleomycin_multicomponent_product"
  sm <- data.frame(source_sampleid=as.character(samples$sampleid),source_sample_name=as.character(samples$sample_name),stringsAsFactors=FALSE)
  sm$named_sampleid <- replace_values(sm$source_sample_name,sm$source_sample_name,policy$sample_display_overrides)
  # Same public-name sanitization as the NCI-60 pipeline; SR needs the panel/cell suffix.
  sm$nci60_2026_sampleid <- paste0("NCI60_",gsub("^_+|_+$","",gsub("[^A-Za-z0-9]+","_",sm$source_sample_name)))
  sm$nci60_2026_sampleid[sm$source_sample_name=="SR"] <- "NCI60_SR_7_19"
  validate_name_maps(drugs,sm)
  list(drugs=drugs,samples=sm)
}
