# Install only into the project's Pixi library, at the tested locked revisions.
.libPaths(.Library)
if (!grepl("/\\.pixi/envs/", normalizePath(.Library)))
  stop("Run pixi run setup; refusing to modify a global R library.")
cfg <- yaml::read_yaml("config/pipeline.yaml")
lock <- jsonlite::fromJSON("renv.lock", simplifyVector=FALSE)
if (as.character(getRversion()) != lock$R$Version) stop("Use the R version in pixi.lock")
Sys.setenv(RENV_CONFIG_CACHE_ENABLED="FALSE", RENV_PATHS_ROOT=file.path(getwd(), ".pixi", "renv"))
options(repos=c(CRAN="https://cloud.r-project.org"), timeout=600)
# Bootstrap a fixed renv tarball; the Conda environment supplies digest already.
if (!requireNamespace("renv", quietly=TRUE) || as.character(packageVersion("renv")) != "1.1.5") {
  archive <- tempfile(fileext=".tar.gz")
  download.file("https://cran.r-project.org/src/contrib/Archive/renv/renv_1.1.5.tar.gz", archive, mode="wb")
  stopifnot(digest::digest(file=archive,algo="sha256") == "3ee996c65989bc5ea3914fbceb3794b00c7daab7eda49dce64048012ea0406b2")
  install.packages(archive, repos=NULL, type="source", lib=.Library)
  unlink(archive)
}
needs_restore <- vapply(lock$Packages, function(record) {
  desc <- tryCatch(utils::packageDescription(record$Package,lib.loc=.Library), error=function(e) NULL)
  !is.list(desc) || !identical(desc$Version,record$Version) ||
    (!is.null(record$RemoteSha) && !identical(desc$RemoteSha,record$RemoteSha))
}, logical(1))
if (any(needs_restore)) renv::restore(project=getwd(),lockfile="renv.lock",library=.Library,
  packages=names(lock$Packages)[needs_restore],prompt=FALSE,transactional=FALSE)
# Validate every locked dependency, not just the two top-level GitHub packages.
for (record in lock$Packages) {
  desc <- utils::packageDescription(record$Package,lib.loc=.Library)
  stopifnot(identical(desc$Version,record$Version))
  if (!is.null(record$RemoteSha)) stopifnot(identical(desc$RemoteSha,record$RemoteSha))
}
manifest <- lapply(names(cfg$packages), function(package) {
  url <- cfg$packages[[package]]
  sha <- sub(".*/commit/", "", url)
  desc <- utils::packageDescription(package,lib.loc=.Library)
  if (!grepl("^[a-f0-9]{40}$",sha) || !identical(desc$RemoteSha,sha)) stop("Configured package revision mismatch: ",package)
  data.frame(package=package,version=desc$Version,source_url=url,resolved_sha=sha)
})
dir.create("data/procdata/metadata",recursive=TRUE,showWarnings=FALSE)
write.table(do.call(rbind,manifest),"data/procdata/metadata/github_packages.tsv",sep="\t",row.names=FALSE,quote=FALSE)
write.table(as.data.frame(installed.packages(lib.loc=.Library)[,c("Package","Version","Built")]),
  "data/procdata/metadata/r_packages.tsv",sep="\t",row.names=FALSE,quote=FALSE)
for (test in c(if(cfg$dataset$id=="HMCL") "test_source_roundtrip.R" else "test_observed_tre.R",
               "test_molecular_annotations.R","test_portable_annotations.R","test_environment.R"))
  local(source(file.path("workflow/tests",test),local=TRUE))
