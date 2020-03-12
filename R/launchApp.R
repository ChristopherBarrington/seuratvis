#' launches the shinyAppDemo app
#'
#' @return shiny application object
#'
#' @example \dontrun{launchApp()}
#'
#' @import shiny
#'
#' @export
launchApp <- function() {
  shinyApp(ui=shinyAppUI, server=shinyAppServer)
}
