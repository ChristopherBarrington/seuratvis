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
  sprintf(fmt='### %s-features_per_cell_filter.ui', id) %>% message()

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

  # record the server(s) to call
  module_environments$features_per_cell_filters$id %<>% c(id) # keep track so other modules can update
  get0(env=module_servers_to_call, x=id) %>% append(sprintf(fmt='%s.server', module)) %>% assign(env=module_servers_to_call, x=id)

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
  session$ns('') %>% sprintf(fmt='### %sfeatures_per_cell_filter.server') %>% message()

  # get environments containing variables to run/configure this object
  collect_environments(id=parent.frame()$id, module='features_per_cell_filter') # provides `seuratvis_env`, `server_env` and `module_env`
  session_server <- get(x='session', env=server_env)

  # react to the minimum input element
  observeEvent(eventExpr=input$min_features, handlerExpr={
    message('### features_per_cell_filter.server-observeEvent-input$min_features')

    # update the reactive
    filtering_parameters.reactions$features_per_cell_min <- round(input$min_features, digits=0)})

  # react to the maximum input element
  observeEvent(eventExpr=input$max_features, handlerExpr={
    message('### features_per_cell_filter.server-observeEvent-input$max_features')

    # update the reactive
    filtering_parameters.reactions$features_per_cell_max <- round(input$max_features, digits=0)})

  # react to the initialisation of the reference min value
  observeEvent(eventExpr=seurat_object.reactions$n_features_values_min, handlerExpr={
    message('### features_per_cell_filter.server-observeEvent-seurat_object.reactions$n_features_values_min')

    # create variables for shorthand
    value <- seurat_object.reactions$n_features_values_min

    # update the ui element(s)
    updateNumericInput(session=session, inputId='min_features', value=value, min=value)})

  # react to the initialisation of the reference max value
  observeEvent(eventExpr=seurat_object.reactions$n_features_values_max, handlerExpr={
    message('### features_per_cell_filter.server-observeEvent-seurat_object.reactions$n_features_values_max')

    # create variables for shorthand
    value <- seurat_object.reactions$n_features_values_max

    # update the ui element(s)
    updateNumericInput(session=session, inputId='max_features', value=value, max=value)})
}
