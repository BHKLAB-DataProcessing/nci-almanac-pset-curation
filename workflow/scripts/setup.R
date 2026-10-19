# Run explicitly with `pixi run setup`; never called by the skeleton workflow.
cfg <- yaml::read_yaml("config/pipeline.yaml")
options(repos = BiocManager::repositories())
manifest <- lapply(cfg$software$github, function(pkg) {
  remotes::install_github(
    repo = pkg$repo, ref = pkg$ref, lib = .Library,
    dependencies = c("Depends", "Imports", "LinkingTo"),
    upgrade = "never", build_vignettes = FALSE, force = TRUE
  )
  desc <- utils::packageDescription(pkg$package, lib.loc = .Library)
  data.frame(package = pkg$package, version = desc$Version,
             repository = pkg$repo, requested_ref = pkg$ref,
             resolved_sha = desc$RemoteSha, stringsAsFactors = FALSE)
})
dir.create("data/procdata/metadata", recursive = TRUE, showWarnings = FALSE)
utils::write.table(do.call(rbind, manifest),
                   "data/procdata/metadata/github_packages.tsv",
                   sep = "\t", row.names = FALSE, quote = FALSE)
