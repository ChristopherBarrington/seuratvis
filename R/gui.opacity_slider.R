#' Set an opacity
#' 
#' Provides a method to set an opacity reactive value
#' 
#' @param id unique name of the element
#' @param label text label of the element
#' 
#' @examples
#' 
#' \dontrun{
#' ui <- fluidPage(opacity_slider.ui(id='page_name'))
#' server <- function(input, output, session) {
#'   callModule(opacity_slider.server, id='page_name')
#' }}
#' 
#' @rdname opacity_slider
#' 
opacity_slider.ui <- function(id, label='Opacity') {
  sprintf(fmt='### %s-opacity_slider.ui', id) %>% message()

  module <- 'opacity_slider'

  # make unique id for this object
  ns <- NS(namespace=id)
  module_ns <- ns(id=module)

  # create an environment in the seuratvis namespace
  e <- new.env()
  e$id <- id
  assign(x=module_ns, val=e, envir=module_environments)

  # record the server(s) to call
  get0(env=module_servers_to_call, x=id) %>% append(sprintf(fmt='%s.server', module)) %>% assign(env=module_servers_to_call, x=id)

  # make ui elements
  sliderInput(inputId=ns(id='opacity_slider'), label=label, min=0.1, max=1, step=0.1, value=1) -> slider

  # return ui element(s)
  tagList(slider)
}

#' React to opacity selection
#' 
#' @details
#' Updates the Seurat object reactive values.
#' 
#' @rdname opacity_slider
#' 
opacity_slider.server <- function(input, output, session) {
  session$ns('') %>% sprintf(fmt='### %sopacity_slider.server') %>% message()
}
