suppressPackageStartupMessages({library(CoreGx);library(data.table)})
source('workflow/scripts/construct_observed_tre.R')
x <- data.table(drug=c('a','a','b','b','c'), dose=c(1,1,2,3,NA), sample=c('x','x','y','z','x'), well=1:5, viab=c(.1,NA,.3,.4,.5))
r<-c('drug','dose'); cc<-c('sample','well'); m<-'viab'
dm<-TREDataMapper(x,rowDataMap=list(r,character()),colDataMap=list(cc,character()),assayMap=list(raw=list(c(r,cc),m)),metadataMap=list())
a<-metaConstruct(dm); b<-construct_observed_tre(x,r,cc,m)
stopifnot(isTRUE(all.equal(a,b)))
z<-x[,.(avg=mean(viab)),by=.(drug,dose,sample)]
a[['summary']]<-z; b[['summary']]<-z
stopifnot(isTRUE(all.equal(a,b)), validObject(b), identical(a[['raw']],b[['raw']]), identical(a[['summary']],b[['summary']]))
cat('Observed-index constructor exactly matches standard constructor and added assay on sparse, NA and repeated-dose fixture.\n')
