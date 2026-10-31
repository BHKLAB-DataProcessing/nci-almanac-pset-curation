# The NCI file mixes tabs and spaces. Parse only the leading NSC token,
# retaining the complete remainder (spaces/punctuation are part of the name).
parse_compound_names <- function(path) {
  lines <- readLines(path,warn=FALSE)
  lines <- lines[nzchar(trimws(lines))]
  good <- grepl("^\\s*[0-9]+\\s+\\S.*$",lines,perl=TRUE)
  if(!all(good)) stop("Unparsed compound-name line(s): ",paste(which(!good),collapse=", "))
  data.table::data.table(
    treatmentid=sub("^\\s*([0-9]+)\\s+.*$","\\1",lines,perl=TRUE),
    treatment_name=trimws(sub("^\\s*[0-9]+\\s+","",lines,perl=TRUE)))
}
