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
  message('### total_umi_per_cell_filter.ui')

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

  module_environments$total_umi_per_cell_filters$ns %<>% c(module_ns)
  module_environments$total_umi_per_cell_filters$id %<>% c(id)

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
  message('### total_umi_per_cell_filter.server')

  # get environments containing variables to run/configure this object
  collect_environments(id=parent.frame()$id, module='total_umi_per_cell_filter') # provides `seuratvis_env`, `server_env` and `module_env`
  session_server <- get(x='session', env=server_env)

  # react to the min_umis input element
  observeEvent(eventExpr=input$min_umis, handlerExpr={
    message('### total_umi_per_cell_filter.server-observeEvent-input$min_umis')

    # update the reactive
    filtering_parameters.reactions$total_umi_per_cell_min <- round(input$min_umis, digits=0)
    seurat_object.reactions$total_umi_per_cell_min <- round(input$min_umis, digits=0)})

  # react to the max_umis input element
  observeEvent(eventExpr=input$max_umis, handlerExpr={
    message('### total_umi_per_cell_filter.server-observeEvent-input$max_umis')

    # update the reactive
    filtering_parameters.reactions$total_umi_per_cell_max <- round(input$max_umis, digits=0)
    seurat_object.reactions$total_umi_per_cell_max <- round(input$max_umis, digits=0)})

  # update UI when Seurat object is loaded
  observeEvent(eventExpr=seurat_object.reactions$seurat, handlerExpr={
    message('### total_umi_per_cell_filter.server-observeEvent-seurat_object.reactions$seurat')

    # create varaibles for shorthand
    seurat <- seurat_object.reactions$seurat
    low <- min(seurat@meta.data$nCount_RNA)
    high <- max(seurat@meta.data$nCount_RNA)

    # update the ui element(s)
    updateNumericInput(session=session, inputId='min_umis', value=low, min=floor(low/50)*50, max=ceiling(high/50)*50)
    updateNumericInput(session=session, inputId='max_umis', value=high, min=floor(low/50)*50, max=ceiling(high/50)*50)

    # update the reactive
    seurat_object.reactions$total_umi_per_cell_min <- low
    seurat_object.reactions$total_umi_per_cell_max <- high})
}
