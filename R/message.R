#' Handle messages to the console
#' 
message <- function(...) {
  if(getOption('seuratvis.verbose', FALSE))
    do.call(base::message, list(...))
  return()
}
