suppressPackageStartupMessages({library(PharmacoGx);library(data.table);library(MultiAssayExperiment)})
setDTthreads(1)
source("workflow/scripts/construct_observed_tre.R")
source("workflow/scripts/named_pset_helpers.R")
source("workflow/scripts/compatibility_v2.R")
source("workflow/scripts/ALMANAC_biomarker_helpers.R")
raw <- CJ(sampleid=c("5|1","1|4"),treatment1dose=c(1e-8,1e-7,1e-6),tech_rep=c("A","B"))
raw[,`:=`(treatment1id="753082",treatment2id="761431",treatment2dose=1e-6,dose_unit="M",NSC1="753082",NSC2="761431",plateid="p",wellid=as.character(.I),viability=seq_len(.N)/10)]
t <- construct_observed_tre(raw,c("treatment1id","treatment2id","treatment1dose","treatment2dose","dose_unit","NSC1","NSC2"),c("sampleid","tech_rep","plateid","wellid"),"viability")
profiles <- unique(raw[,.(sampleid,tech_rep,treatment1id,treatment2id,NSC1,NSC2)])
profiles[,mean_Bliss_observed_excess:=seq_len(.N)/100]
t[["profiles"]] <- profiles
sd <- data.frame(sampleid=c("5|1","1|4"),sample_name=c("MCF7","A549/ATCC"),row.names=c("5|1","1|4"))
dd <- data.frame(treatmentid=c("753082","761431"),treatment_name="Vemurafenib",row.names=c("753082","761431"))
se <- SummarizedExperiment(assays=list(value=matrix(1:6,3,2,dimnames=list(letters[1:3],sd$sampleid))),colData=S4Vectors::DataFrame(sampleid=sd$sampleid,batchid="fixture",row.names=sd$sampleid))
metadata(se)$annotation <- "rnaseq.composite"
mae <- MultiAssayExperiment(experiments=list(rnaseq.composite=se),colData=S4Vectors::DataFrame(sd))
p <- PharmacoSet2(name="fixture",sample=sd,treatment=dd,molecularProfiles=mae,treatmentResponse=t,curation=list(sample=data.frame(sampleid=sd$sampleid,sampleid_dataset=sd$sampleid,row.names=sd$sampleid),treatment=data.frame(treatmentid=dd$treatmentid,treatmentid_dataset=dd$treatmentid,row.names=dd$treatmentid),tissue=data.frame()),datasetType="sensitivity")
d <- data.frame(nsc=dd$treatmentid,named_treatmentid=paste0("Vemurafenib [NSC ",dd$treatmentid,"]"),preferred_name="Vemurafenib",release_cas=NA_character_,identity_review_flags=NA_character_)
s <- data.frame(source_sampleid=sd$sampleid,named_sampleid=c("MCF-7","A-549"),nci60_2026_sampleid=c("NCI60_MCF7","NCI60_A549_ATCC"))
v1 <- name_almanac_pset(p,d,s)
before <- compatibility_scientific_fingerprint(v1)
v2 <- almanac_compatibility_v2(v1)
stopifnot(identical(before,compatibility_scientific_fingerprint(v2)),nrow(v2@treatment)==1L,
  identical(drugNames(v2),"Vemurafenib"),v2@treatment$treatmentid_dataset=="753082///761431",
  nrow(v2@treatmentResponse[["raw"]])==12L,nrow(v2@treatmentResponse[["profiles"]])==4L,
  identical(v1@treatmentResponse@rowData$treatment1id,rep("Vemurafenib [NSC 753082]",nrow(v1@treatmentResponse@rowData))),
  identical(assayNames(v2@molecularProfiles[["rnaseq.comp"]]),"exprs"),!"rna" %in% names(v2@molecularProfiles))
stopifnot(inherits(try(almanac_pair_profiles(v2,c("Vemurafenib","Vemurafenib")),silent=TRUE),"try-error"))
selected <- almanac_pair_profiles(v2,c("Vemurafenib","Vemurafenib"),c("753082","761431"))
stopifnot(nrow(selected)==4L,all(selected$NSC1=="753082"),all(selected$NSC2=="761431"))
stopifnot(inherits(try(almanac_pair_profiles(v2,c("Vemurafenib","Vemurafenib"),c("1","761431")),silent=TRUE),"try-error"))
f <- tempfile(fileext=".rds");saveRDS(v2,f);z<-readRDS(f)
stopifnot(identical(before,compatibility_scientific_fingerprint(z)));unlink(f)
cat("PASS: canonical treatment grouping retains all 12 wells/4 profiles and two source NSCs; exact assays/indices/molecular values, input-object preservation, RNA labels, ambiguity guard, source selection and save/reload\n")
