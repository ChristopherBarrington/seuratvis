#' Select element of vecto by preference
#' 
#' Returns the element of a vector using a ranked preference list
#' 
#' @param x vector from which to choose
#' @param preferences ranked vector of terms to serach for in \code{x}
#' @param default if no \code{preferences} are found, the \code{default} element of \code{x} is returned
#' 
#' @return an element of \code{x}
#' 
preferred_choice <- function(x, preferences, default=1) {
  preferences %<>% append(pluck(x, default))  
  x %>% str_c(collapse='|') %>% str_subset(string=preferences) %>% head(n=1)
}
