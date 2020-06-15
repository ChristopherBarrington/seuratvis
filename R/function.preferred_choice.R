#' Select element of vecto by preference
#' 
#' Returns the element of a vector using a ranked preference list
#' 
#' @param x vector from which to choose
#' @param preferences ranked vector of terms to serach for in \code{x}
#' @param default if no \code{preferences} are found, the \code{default}-th element of \code{x} is returned
#' 
#' @return an element of \code{x}
#' 
preferred_choice <- function(x, preferences, default=1) {
  if(length(x)==0)
    stop('no elements from which to choose!')
  
  # if no `preferences` are in `x`, use the first element of `x`
  preferences %<>% append(pluck(x, default))

  # return the first match in ranked `preferences` that is in `x`
  x %>%
    str_c(collapse='|') %>%
    sprintf(fmt='^(%s)$') %>%
    str_subset(string=preferences) %>%
    head(n=1)
}
