# packageDescription carries an installation-specific `file` attribute and may
# carry LibPath. Keep version/RemoteSha provenance while removing local paths.
portable_annotations <- function(x) {
  if (inherits(x, "sessionInfo")) {
    for (name in intersect(c("BLAS", "LAPACK"), names(x))) {
      if (is.character(x[[name]])) x[[name]] <- basename(x[[name]])
    }
  }
  if (inherits(x, "packageDescription")) {
    attr(x, "file") <- NULL
    x$LibPath <- NULL
  }
  if (is.list(x)) {
    for (i in seq_along(x)) x[i] <- list(portable_annotations(x[[i]]))
  }
  x
}
