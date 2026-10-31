# CoreGx 2.7.3 LongTable() expands all row/column combinations before
# discarding unobserved cells. Construct the same raw-only index by joins.
# This adapter uses the installed CoreGx class layout; equivalence to the
# public constructor is checked on fixtures, including repeated dose wells.
construct_observed_tre <- function(raw, row_keys, col_keys, measures) {
  ids <- c(row_keys, col_keys)
  stopifnot(!anyDuplicated(raw, by=ids))
  rd <- unique(raw[, row_keys, with=FALSE]); setorderv(rd, row_keys)
  cd <- unique(raw[, col_keys, with=FALSE]); setorderv(cd, col_keys)
  rd[, rowKey := .I]; cd[, colKey := .I]
  values <- copy(raw[, c(ids, measures), with=FALSE])
  setkeyv(values, ids); values[, .raw := .I]
  index <- values[, c(ids, ".raw"), with=FALSE]
  index <- merge(index, rd, by=row_keys, sort=FALSE)
  index <- merge(index, cd, by=col_keys, sort=FALSE)
  stopifnot(nrow(index) == nrow(raw), !anyNA(index$rowKey), !anyNA(index$colKey))
  index <- index[, .(colKey, rowKey, .raw)]; setkeyv(index, ".raw")
  values[, (ids) := NULL]; setkeyv(values, ".raw")
  intern <- CoreGx:::immutable(list(rowIDs=row_keys, rowMeta=character(),
    colIDs=col_keys, colMeta=character(), assayKeys=list(raw=ids), assayIndex=index))
  rd[, .rownames := do.call(paste, c(.SD, sep=":")), .SDcols=row_keys]
  cd[, .colnames := do.call(paste, c(.SD, sep=":")), .SDcols=col_keys]
  setkeyv(rd, "rowKey"); setkeyv(cd, "colKey")
  lt <- methods::new("LongTable", rowData=rd, colData=cd,
    assays=list(raw=values), metadata=list(), .intern=intern)
  tre <- CoreGx::TreatmentResponseExperiment(lt)
  methods::validObject(tre)
  tre
}
