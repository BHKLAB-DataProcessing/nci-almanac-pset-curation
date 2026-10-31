source('workflow/scripts/check_environment.R')
packages<-yaml::read_yaml('config/pipeline.yaml')$packages
stopifnot(check_curation_environment(packages))
stopifnot(inherits(try(check_curation_environment(packages,library='/global/R/library'),silent=TRUE),'try-error'))
manifest<-data.table::fread('data/procdata/metadata/github_packages.tsv')
manifest[package=='PharmacoGx', resolved_sha:='unapproved-revision']
path<-tempfile(fileext='.tsv');data.table::fwrite(manifest,path,sep='\t')
stopifnot(inherits(try(check_curation_environment(packages,manifest_path=path),silent=TRUE),'try-error'))
unlink(path)
cat('PASS: recorded Pixi environment accepted; global library and mismatched framework revision rejected.\n')
