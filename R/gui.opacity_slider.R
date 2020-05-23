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
  message('### opacity_slider.ui')

  module <- 'opacity_slider'

  # make unique id for this object
  ns <- NS(namespace=id)
  module_ns <- ns(id=module)

  # create an environment in the seuratvis namespace
  e <- new.env()
  e$id <- id
  assign(x=module_ns, val=e, envir=module_environments)

  module_environments$opacity_sliders$ns %<>% c(module_ns)
  module_environments$opacity_sliders$id %<>% c(id)

  # make ui elements
  ## if a label switch is required, make one
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
  message('### opacity_slider.server')

  # get environments containing variables to run/configure this object
  collect_environments(id=parent.frame()$id, module='opacity_slider') # provides `seuratvis_env`, `server_env` and `module_env`
  session_server <- get(x='session', env=server_env)

  # react to the opacity slider
  observeEvent(eventExpr=input$opacity_slider, handlerExpr={
    message('### opacity_slider.server-observeEvent-input$opacity_slider')

    # update the reactive
    seurat_object.reactions$opacity <- input$opacity_slider})

  # update UI when Seurat object is loaded
  observeEvent(eventExpr=seurat_object.reactions$seurat, handlerExpr={
    message('### opacity_slider.server-observeEvent-seurat_object.reactions$seurat')

    # update the reactive
    seurat_object.reactions$opacity <- 1})
}
