#' Specify total UMIs per cell limits
#' 
#' Provides a way to define limits on total UMI per cell
#' 
#' @param id unique name of the element
#' @param label text label of the element
#' @param low,high include inputs for lower and upper thresholds, defaults to min/max of Seurat object
#' 
#' @examples
#' 
#' \dontrun{
#' ui <- fluidPage(total_umi_per_cell_filter.ui(id='page_name'))
#' server <- function(input, output, session) {
#'   callModule(total_umi_per_cell_filter.server, id='page_name')
#' }} 
#' 
#' @rdname total_umi_per_cell_filter
#' 
total_umi_per_cell_filter.ui <- function(id, label='UMIs per cell', low=TRUE, high=TRUE) {
  sprintf(fmt='### %s-total_umi_per_cell_filter.ui', id) %>% message()

  module <- 'total_umi_per_cell_filter'

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
  module_environments$total_umi_per_cell_filters$id %<>% c(id) # keep track so other modules can update
  get0(env=module_servers_to_call, x=id) %>% append(sprintf(fmt='%s.server', module)) %>% assign(env=module_servers_to_call, x=id)

  # make ui elements
  low_ui <- high_ui <- NULL

  ## if a minium box is required, make one
  if(low)
    div(tags$h6('Minimum', style='display: inline;'),
        numericInput(inputId=ns('min_umis'), label=NULL, value=0, step=50, width='100%')) -> low_ui

  ## if a max_umis box is required, make one
  if(high)
    div(tags$h6('Maximum', style='display: inline;'),
        numericInput(inputId=ns('max_umis'), label=NULL, value=0, step=50, width='100%')) -> high_ui

  # return ui element(s)
  tagList(tags$label(label), br(), low_ui, high_ui)
}

#' React to total UMI per cell thresholds
#' 
#' @details
#' Updates the Seurat object reactive values.
#' 
#' @rdname total_umi_per_cell_filter
#'
total_umi_per_cell_filter.server <- function(input, output, session) {
  session$ns('') %>% sprintf(fmt='### %stotal_umi_per_cell_filter.server') %>% message()

  # get environments containing variables to run/configure this object
  collect_environments(id=parent.frame()$id, module='total_umi_per_cell_filter') # provides `seuratvis_env`, `server_env` and `module_env`
  session_server <- get(x='session', env=server_env)

  # react to the min_umis input element
  observeEvent(eventExpr=input$min_umis, handlerExpr={
    # send a message
    session$ns('') %>% sprintf(fmt='### %stotal_umi_per_cell_filter.server-observeEvent-input$min_umis [%s]', input$min_umis) %>% message('')

    # update the reactive
    filtering_parameters.reactions$total_umi_per_cell_min <- round(input$min_umis, digits=0)})

  # react to the max_umis input element
  observeEvent(eventExpr=input$max_umis, handlerExpr={
    # send a message
    session$ns('') %>% sprintf(fmt='### %stotal_umi_per_cell_filter.server-observeEvent-input$max_umis [%s]', input$max_umis) %>% message('')

    # update the reactive
    filtering_parameters.reactions$total_umi_per_cell_max <- round(input$max_umis, digits=0)})

  # react to the initialisation of the reference min value
  observeEvent(eventExpr=seurat_object.reactions$n_umi_values_min, handlerExpr={
    # send a message
    session$ns('') %>% sprintf(fmt='### %s-observeEvent-seurat_object.reactions$n_umi_values_min [%s]', seurat_object.reactions$n_umi_values_min) %>% message('')

    # create variables for shorthand
    value <- seurat_object.reactions$n_umi_values_min

    # update the ui element(s)
    updateNumericInput(session=session, inputId='min_umis', value=value, min=value)})

  # react to the initialisation of the reference max value
  observeEvent(eventExpr=seurat_object.reactions$n_umi_values_max, handlerExpr={
    # send a message
    session$ns('') %>% sprintf(fmt='### %s-observeEvent-seurat_object.reactions$n_umi_values_max [%s]', seurat_object.reactions$n_umi_values_max) %>% message('')

    # create variables for shorthand
    value <- seurat_object.reactions$n_umi_values_max

    # update the ui element(s)
    updateNumericInput(session=session, inputId='max_umis', value=value, max=value)})
}
