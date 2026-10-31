suppressPackageStartupMessages({library(CoreGx);library(data.table)})
source('workflow/qc/source_roundtrip.R')
x<-CJ(sampleid=c('s1','s2'),treatment1dose=c(.1,1,10),tech_rep=c('r1','r2'))
x[,`:=`(treatment1id='drug1',treatment2id=NA_character_,treatment2dose=NA_real_,dose_unit='uM',plateid='plate',viability=seq_len(.N)/20)]
x[3,viability:=NA_real_]
r<-c('treatment1id','treatment2id','treatment1dose','treatment2dose','dose_unit');cc<-c('sampleid','tech_rep','plateid')
dm<-TREDataMapper(x,rowDataMap=list(r,character()),colDataMap=list(cc,character()),assayMap=list(raw=list(c(r,cc),'viability')),metadataMap=list())
t<-metaConstruct(dm);p<-tempfile(fileext='.rds');saveRDS(t,p);t<-readRDS(p);unlink(p)
cfg<-list(dataset=list(id='HMCL'))
y<-t[['raw']]
stopifnot(nrow(y)==12L,all(validate_source_roundtrip(y,cfg,src=x,output=NULL)$pass))
stopifnot(all(validate_source_roundtrip(y[12:1],cfg,src=x,output=NULL)$pass))
bad<-copy(y);bad[1,viability:=.999]
stopifnot(inherits(try(validate_source_roundtrip(bad,cfg,src=x,output=NULL),silent=TRUE),'try-error'))
stopifnot(inherits(try(validate_source_roundtrip(y[c(1,7)],cfg,src=x,output=NULL),silent=TRUE),'try-error'))
cat('PASS: serialized 2 samples x 3 doses x 2 replicates retain every value; row reordering allowed; collapse and value corruption rejected.\n')
