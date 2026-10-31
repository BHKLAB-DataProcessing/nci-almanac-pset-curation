# An API is mutable. Live mode permits new responses; replay verifies exact bytes.
annotationdb_read_response <- function(url,cache,key,cfg) {
  mode <- cfg$annotationdb$mode
  if (is.null(mode) || !mode %in% c("live","replay")) stop("annotationdb.mode must be live or replay")
  expected <- NULL
  if(mode=="replay") {
    pins <- jsonlite::fromJSON("config/annotationdb-checksums.json",simplifyVector=FALSE)
    expected <- pins[[key]]
    if(is.null(expected)) stop("No pinned AnnotationDB response for this query; choose an explicit live refresh and review the result: ",key)
  }
  refresh <- isTRUE(cfg$annotationdb$refresh)
  if(mode=="replay" && refresh) stop("Replay and refresh are mutually exclusive")
  if(!file.exists(cache) || refresh) {
    response <- httr2::request(url) |> httr2::req_timeout(180) |> httr2::req_retry(max_tries=3) |> httr2::req_perform()
    bytes <- httr2::resp_body_raw(response)
    actual <- digest::digest(bytes,algo="sha256",serialize=FALSE)
    if(!is.null(expected) && !identical(actual,expected)) stop("AnnotationDB changed; archived response required for exact replay: ",key)
    # Validate JSON before atomically promoting a response.
    jsonlite::fromJSON(rawToChar(bytes),simplifyVector=FALSE)
    temp <- paste0(cache,".part");writeBin(bytes,temp)
    if(!file.rename(temp,cache)) stop("Cannot save AnnotationDB response")
  }
  actual <- digest::digest(file=cache,algo="sha256")
  if(!is.null(expected) && !identical(actual,expected)) stop("Cached AnnotationDB response differs from replay checksum: ",key)
  data.frame(key=key,sha256=actual,mode=mode,stringsAsFactors=FALSE)
}
