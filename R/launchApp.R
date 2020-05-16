#' launches the shinyAppDemo app
#'
#' @return shiny application object
#'
#' @import shiny
#'
#' @export
launchApp <- function() {
  shinyApp(ui=shinyAppUI, server=shinyAppServer, options=list(quiet=!Sys.info()['user'] %in% c('chris')))
}
