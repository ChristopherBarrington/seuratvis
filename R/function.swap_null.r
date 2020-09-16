#'
#'
swap_null <- function(x, default)
  ifelse(is.null(x), default, x)
