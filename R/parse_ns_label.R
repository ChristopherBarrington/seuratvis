#' Parse NS label
#' 
#' Remove the interim namespace label
#' 
#' @details Removes the intermediate namespace tag where the same module is use in multiple parts of the same tab. For example \code{tab-module-id} becomes \code{tab-id}.
#' 
parse_ns_label <- function(x)
  str_replace(string=x, pattern='-.*-', replacement='-')
