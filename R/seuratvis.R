#' launches the seuratvis app
#'
#' @return shiny application object
#'
#' @import shiny
#'
#' @examples
#' 
#' \dontrun{
#' seuratvis()}
#' 
#' @export
#' 
seuratvis <- function()
  shinyApp(ui=shinyAppUI, server=shinyAppServer)
