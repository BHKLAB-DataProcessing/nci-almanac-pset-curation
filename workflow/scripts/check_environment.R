# Keep the global R installation and unpublished source checkouts out of runs.
check_curation_environment <- function(packages,
    manifest_path="data/procdata/metadata/github_packages.tsv", library=.Library) {
  library <- normalizePath(library, mustWork=FALSE)
  if (!grepl("/\\.pixi/envs/", library))
    stop("Use the project's Pixi R environment; default/global R is not supported: ", library)
  if (!file.exists(manifest_path)) stop("Run pixi run setup to record exact GitHub package revisions")
  manifest <- data.table::fread(manifest_path)
  for (package_name in names(packages)) {
    expected <- manifest[manifest$package == package_name]
    description <- utils::packageDescription(package_name)
    package_path <- normalizePath(find.package(package_name), mustWork=TRUE)
    if (nrow(expected) != 1L || !is.list(description) || is.null(description$RemoteSha) ||
        !identical(as.character(description$RemoteSha), as.character(expected$resolved_sha)) ||
        !identical(as.character(packages[[package_name]]), as.character(expected$source_url)) ||
        !startsWith(package_path, paste0(library, "/")))
      stop("Package environment does not match configured, recorded GitHub revision: ", package_name,
           ". Run pixi run setup in an idle environment before continuing.")
  }
  invisible(TRUE)
}
