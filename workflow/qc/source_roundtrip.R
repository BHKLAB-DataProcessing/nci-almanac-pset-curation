# Validate stored source measurements through a public raw-assay accessor.
# Counts and S4 validity alone cannot rule out dose/replicate index loss.
validate_source_roundtrip <- function(raw, cfg, src=NULL, output='data/procdata/metadata/source_roundtrip.tsv') {
  if (is.null(src)) src <- readRDS('data/procdata/response.rds')$response
  # Named public identifiers are labels; compare immutable source NSCs/context IDs.
  if(cfg$dataset$id=='NCI_ALMANAC') {
    raw <- copy(raw)
    raw[, `:=`(treatment1id=as.character(NSC1),treatment2id=as.character(NSC2))]
    raw[treatment2id=="", treatment2id:=NA_character_]
    source('workflow/scripts/named_pset_helpers.R',local=TRUE)
    source('workflow/scripts/almanac_source_names.R',local=TRUE)
    b <- readRDS('data/procdata/response.rds')
    m <- almanac_source_name_maps(b$sample,b$treatment,
      names_file='data/rawdata/downloads/combocompoundnames_small/payload',
      structures_file='data/rawdata/downloads/combocompoundset/payload')
    raw[, sampleid := map_ids(sampleid,setNames(m$samples$source_sampleid,m$samples$named_sampleid))]
  }
  row_keys <- c('treatment1id','treatment2id','treatment1dose','treatment2dose','dose_unit')
  col_keys <- c('sampleid','tech_rep','plateid')
  if(cfg$dataset$id=='GDSC_MATRIX')col_keys<-c(col_keys,'matrixid')
  if(cfg$dataset$id=='GDSC_ANCHOR')row_keys<-c(row_keys,'library_doses','anchor_doses')
  if(cfg$dataset$id=='NCI_ALMANAC')row_keys<-c(row_keys,'NSC1','NSC2')
  if(cfg$dataset$id!='HMCL')col_keys<-c(col_keys,'wellid')
  keys<-c(row_keys,col_keys)
  measures<-setdiff(names(src)[vapply(src,function(x)is.numeric(x)||is.logical(x),logical(1))],c(keys,'dose_index'))
  required<-c(keys,measures)
  stopifnot(all(required%in%names(raw)),all(keys%in%names(src)))
  checks<-data.table(check=c('row_count','source_keys_unique','object_keys_unique'),
    pass=c(nrow(src)==nrow(raw),!anyDuplicated(src,by=keys),!anyDuplicated(raw,by=keys)))
  if(!all(checks$pass)) {
    if(!is.null(output))fwrite(checks,output,sep='\t')
    stop('Source/assay row count or physical-context key uniqueness failed')
  }
  si<-copy(src[,..keys]);ri<-copy(raw[,..keys])
  si[,validation_order:=.I];ri[,validation_order:=.I]
  setorderv(si,keys);setorderv(ri,keys)
  so<-si$validation_order;ro<-ri$validation_order
  rm(si,ri)
  for(n in required) {
    a<-src[[n]][so];b<-raw[[n]][ro]
    if((is.numeric(a)||is.logical(a))&&(is.numeric(b)||is.logical(b))) {
      a<-as.numeric(a);b<-as.numeric(b)
    } else {a<-as.character(a);b<-as.character(b)}
    checks<-rbind(checks,data.table(check=paste0(if(n%in%keys)'key:' else 'measurement:',n),pass=identical(a,b)))
  }
  checks[,`:=`(source_rows=nrow(src),object_rows=nrow(raw),dataset=cfg$dataset$id,
    CoreGx_version=as.character(utils::packageVersion('CoreGx')),
    CoreGx_sha=utils::packageDescription('CoreGx')$RemoteSha)]
  if(!is.null(output))fwrite(checks,output,sep='\t')
  if(!all(checks$pass))stop('Source-to-object key/value parity failed: ',paste(checks[pass==FALSE,check],collapse=', '))
  message('PASS source round trip: ',nrow(src),' rows, ',length(keys),' dose/replicate keys, ',length(measures),' source measurement columns.')
  checks
}
