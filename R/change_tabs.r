#'
#'
go_to_config <- function(inputId='left_sidebar', session)
  go_to_tab(tab='configuration_tab', session=session, inputId=inputId)

#'
#' 
go_to_tab <- function(tab, inputId, session)
  updateNavbarPage(session=session, inputId=inputId, selected=tab)
