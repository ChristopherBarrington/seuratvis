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
  message('### opacity_slider.ui')

  module <- 'point_size_slider'

  # make unique id for this object
  ns <- NS(namespace=id)
  module_ns <- ns(id=module)

  # create an environment in the seuratvis namespace
  e <- new.env()
  e$id <- id
  assign(x=module_ns, val=e, envir=module_environments)

  module_environments$point_size_sliders$ns %<>% c(module_ns)
  module_environments$point_size_sliders$id %<>% c(id)

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
point_size_slider.server <- function(input, output, session) {
  message('### point_size_slider.server')

  # get environments containing variables to run/configure this object
  collect_environments(id=parent.frame()$id, module='point_size_slider') # provides `seuratvis_env`, `server_env` and `module_env`
  session_server <- get(x='session', env=server_env)

  # react to the point size slider
  observeEvent(eventExpr=input$point_size_slider, handlerExpr={
    message('### point_size_slider.server-observeEvent-input$point_size_slider')

    # update the reactive
    seurat_object.reactions$point_size <- input$point_size_slider})

  # update UI when Seurat object is loaded
  observeEvent(eventExpr=seurat_object.reactions$seurat, handlerExpr={
    sprintf(fmt='### point_size_slider.server-observeEvent-seurat_object.reactions$seurat [%s]', seurat_object.reactions$formatted.project.name) %>% message()

    # update the reactive
    seurat_object.reactions$point_size <- 0.6})
}
