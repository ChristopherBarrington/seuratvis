#' Specify unique features per cell limits
#' 
#' Provides a way to define limits on number of features identified in a cell
#' 
#' @param id unique name of the element
#' @param label text label of the element
#' @param low,high include inputs for lower and upper thresholds, defaults to min/max of Seurat object
#' 
#' @examples
#' 
#' \dontrun{
#' ui <- fluidPage(features_per_cell_filter.ui(id='page_name'))
#' server <- function(input, output, session) {
#'   callModule(features_per_cell_filter.server, id='page_name')
#' }}
#' 
#' @rdname features_per_cell_filter
#' 
features_per_cell_filter.ui <- function(id, label='Features per cell', low=TRUE, high=TRUE) {
  message('### features_per_cell_filter.ui')

  module <- 'features_per_cell_filter'

  # make unique id for this object
  ns <- NS(namespace=id)
  module_ns <- ns(id=module)

  # create an environment in the seuratvis namespace
  e <- new.env()
  e$id <- id
  e$low <- low
  e$high <- high
  assign(x=module_ns, val=e, envir=module_environments)

  module_environments$features_per_cell_filters$ns %<>% c(module_ns)
  module_environments$features_per_cell_filters$id %<>% c(id)

  # make ui elements
  low_ui <- high_ui <- NULL

  ## if a minium box is required, make one
  if(low)
    div(tags$h6('Minimum', style='display: inline;'),
        numericInput(inputId=ns('min_features'), label=NULL, value=0, step=50, width='100%')) -> low_ui

  ## if a maximum box is required, make one
  if(high)
    div(tags$h6('Maximum', style='display: inline;'),
        numericInput(inputId=ns('max_features'), label=NULL, , value=0, step=50, width='100%')) -> high_ui

  # return ui element(s)
  tagList(tags$label(label), br(), low_ui, high_ui)
}

#' React to features per cell thresholds
#' 
#' @details
#' Updates the Seurat object reactive values.
#' 
#' @rdname features_per_cell_filter
#'
features_per_cell_filter.server <- function(input, output, session) {
  message('### features_per_cell_filter.server')

  # get environments containing variables to run/configure this object
  collect_environments(id=parent.frame()$id, module='features_per_cell_filter') # provides `seuratvis_env`, `server_env` and `module_env`
  session_server <- get(x='session', env=server_env)

  # react to the minimum input element
  observeEvent(eventExpr=input$min_features, handlerExpr={
    message('### features_per_cell_filter.server-observeEvent-input$min_features')

    # update the reactive
    filtering_parameters.reactions$features_per_cell_min <- round(input$min_features, digits=0)
    seurat_object.reactions$features_per_cell_min <- round(input$min_features, digits=0)})

  # react to the maximum input element
  observeEvent(eventExpr=input$max_features, handlerExpr={
    message('### features_per_cell_filter.server-observeEvent-input$max_features')

    # update the reactive
    filtering_parameters.reactions$features_per_cell_max <- round(input$max_features, digits=0)
    seurat_object.reactions$features_per_cell_max <- round(input$max_features, digits=0)})

  # update UI when Seurat object is loaded
  observeEvent(eventExpr=seurat_object.reactions$seurat, handlerExpr={
    sprintf(fmt='### features_per_cell_filter.server-observeEvent-seurat_object.reactions$seurat [%s]', seurat_object.reactions$formatted.project.name) %>% message()

    # create varaibles for shorthand
    seurat <- seurat_object.reactions$seurat
    low <- min(seurat@meta.data$nFeature_RNA)
    high <- max(seurat@meta.data$nFeature_RNA)

    # update the ui element(s)
    updateNumericInput(session=session, inputId='min_features', value=low, min=floor(low/50)*50, max=ceiling(high/50)*50)
    updateNumericInput(session=session, inputId='max_features', value=high, min=floor(low/50)*50, max=ceiling(high/50)*50)

    # update the reactive
    seurat_object.reactions$features_per_cell_min <- low
    seurat_object.reactions$features_per_cell_max <- high})
}
