#' Set point sizes
#' 
#' Provides a method to set a point size reactive value
#' 
#' @param id unique name of the element
#' @param label text label of the element
#' 
#' @examples
#' 
#' \dontrun{
#' ui <- fluidPage(point_size_slider.ui(id='page_name'))
#' server <- function(input, output, session) {
#'   callModule(point_size_slider.server, id='page_name')
#' }}
#' 
#' @rdname point_size_slider
#' 
point_size_slider.ui <- function(id, label='Point size') {
  sprintf(fmt='### %s-point_size_slider.ui', id) %>% message()

  module <- 'point_size_slider'

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
  ## if a label switch is required, make one
  sliderInput(inputId=ns(id='point_size_slider'), label=label, min=0.1, max=3, step=0.1, value=0.6) -> slider

  # return ui element(s)
  tagList(slider)
}

#' React to point size selection
#' 
#' @details
#' Updates the Seurat object reactive values.
#' 
#' @rdname point_size_slider
#' 
point_size_slider.server <- function(input, output, session, ...) {
  session$ns('') %>% sprintf(fmt='### %spoint_size_slider.server') %>% message()
}
