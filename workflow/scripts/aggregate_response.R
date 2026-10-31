# Preserve singleton dose contexts directly. Aggregate only when needed, using
# data.table's grouped mean after applying the documented finite-value policy.
summarize_response <- function(x, keys, measures) {
  out <- data.table::copy(x[, c(keys, measures), with=FALSE])
  for (column in measures) {
    bad <- which(!is.finite(out[[column]]))
    if (length(bad)) data.table::set(out, i=bad, j=column, value=NA)
  }
  if (anyDuplicated(out, by=keys)) {
    counts <- out[, .(n_observations=.N), by=keys]
    out <- out[, lapply(.SD, mean, na.rm=TRUE), by=keys, .SDcols=measures]
    out <- merge(out, counts, by=keys, sort=FALSE)
    for (column in measures) {
      bad <- which(!is.finite(out[[column]]))
      if (length(bad)) data.table::set(out, i=bad, j=column, value=NA)
    }
  } else out[, n_observations:=1L]
  out
}
