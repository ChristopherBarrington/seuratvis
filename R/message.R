#' Handle messages to the console
#' 
message <- function(...) {
  if(Sys.getenv('USER')=='chris') 
    do.call(base::message, list(...))
  return()
}
