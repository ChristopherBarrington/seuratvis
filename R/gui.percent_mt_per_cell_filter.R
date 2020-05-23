#' Specify maximum mitochondrial proportion
#' 
#' Provides a way to define limit on proportion of mitochondrial UMI per cell
#' 
#' @param id unique name of the element
#' @param label text label of the element
#' @param target_var name of column in \code{@meta.data} that is percentage mitochondria
#' 
#' @examples
#' 
#' \dontrun{
#' ui <- fluidPage(percent_mt_per_cell_filter.ui(id='page_name'))
#' server <- function(input, output, session) {
#'   callModule(percent_mt_per_cell_filter.server, id='page_name')
#' }}
#' 
#' @rdname percent_mt_per_cell_filter
#' 
percent_mt_per_cell_filter.ui <- function(id, label='Proportion Mt', target_var='percent_mt') {
  message('### features_per_cell_filter.ui')

  module <- 'percent_mt_per_cell_filter'

  # make unique id for this object
  ns <- NS(namespace=id)
  module_ns <- ns(id=module)

  # create an environment in the seuratvis namespace
  e <- new.env()
  e$id <- id
  e$target_var <- target_var
  assign(x=module_ns, val=e, envir=module_environments)

  module_environments$percent_mt_per_cell_filters$ns %<>% c(module_ns)
  module_environments$percent_mt_per_cell_filters$id %<>% c(id)

  # make ui elements
  textInput(inputId=ns('max_percent_mt'), label='Maximum', placeholder='Maximum proportion Mt') -> high_ui

  # return ui element(s)
  tagList(tags$label(label), br(), high_ui)
}

#' React to maximum proportion mitochondria per cell threshold
#' 
#' @details
#' Updates the Seurat object reactive values.
#' 
#' @rdname percent_mt_per_cell_filter
#'
percent_mt_per_cell_filter.server <- function(input, output, session) {
  message('### percent_mt_per_cell_filter.server')

  # get environments containing variables to run/configure this object
  collect_environments(id=parent.frame()$id, module='percent_mt_per_cell_filter') # provides `seuratvis_env`, `server_env` and `module_env`
  session_server <- get(x='session', env=server_env)

  # react to the maximum input element
  observeEvent(eventExpr=input$max_percent_mt, handlerExpr={
    message('### percent_mt_per_cell_filter.server-observeEvent-input$max_percent_mt')

    # update the reactive
    seurat_object.reactions$percent_mt_per_cell_max <- input$max_percent_mt %>% as.numeric() %>% add(0.05) %>% round(digits=1)})

  # update UI when Seurat object is loaded
  observeEvent(eventExpr=seurat_object.reactions$seurat, handlerExpr={
    message('### percent_mt_per_cell_filter.server-observeEvent-seurat_object.reactions$seurat')

    # create varaibles for shorthand
    seurat <- seurat_object.reactions$seurat
    high <- max(seurat@meta.data$percent_mt) %>% add(0.05) %>% round(digits=1)

    # update the ui element(s)
    updateTextInput(session=session, inputId='max_percent_mt', placeholder=high)

    # update the reactive
    seurat_object.reactions$percent_mt_per_cell_max <- high
    seurat_object.reactions$percent_mt_target_var <- module_env$target_var
    seurat_object.reactions$percent_mt <- FetchData(object=seurat, vars=module_env$target_var) %>% set_names('percent_mt')})
}
