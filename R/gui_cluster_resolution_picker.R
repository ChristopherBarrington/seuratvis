#' Select cluster resolution
#' 
#' Provides a method to select a cluster resolution to plot
#' 
#' @param id unique name of the element
#' @param label text label of the element
#' 
#' @examples
#' 
#' \dontrun{
#' ui <- fluidPage(cluster_resolution_picker.ui(id='page_name'))
#' server <- function(input, output, session) {
#'   callModule(update_selected_cluster_resolution.server, id='page_name')
#' }}
#' 
#' 
#' @rdname cluster_resolution_picker
#' 
cluster_resolution_picker.ui <- function(id, label='Cluster resolutions') {
  message('### cluster_resolution_picker.ui')

  module <- 'cluster_resolution_picker'

  # make unique id for this object
  ns <- NS(namespace=id, id=module)

  # create an environment in the seuratvis namespace
  e <- new.env()
  e$id <- id
  assign(x=ns, val=e, envir=module_environments)

  module_environments$cluster_resolution_pickers$ns %<>% c(ns)
  module_environments$cluster_resolution_pickers$id %<>% c(id)

  # return ui element(s)
  selectInput(inputId=ns, label=label,
              choices='seurat_clusters', selected='seurat_clusters',
              multiple=FALSE)
}

#' React to a cluster resolution choice
#' 
#' @details
#' Updates the Seurat object reactive values. Could use \code{SetIdent} instead?
#' 
#' @rdname cluster_resolution_picker
#' 
update_selected_cluster_resolution.server <- function(input, output, session, i) {
  message('### update_selected_cluster_resolution.server')

  # get environments containing variables to run/configure this object
  collect_environments(id=parent.frame()$id, module='cluster_resolution_picker') # provides `seuratvis_env`, `server_env` and `module_env`
  session <- get(x='session', env=server_env)

  # react to the cluster resolution selection
  observeEvent(eventExpr=input$cluster_resolution_picker, handlerExpr={
    r <- input$cluster_resolution_picker
    seurat_object.reactions$selected_cluster_resolution <- r
    seurat_object.reactions$selected_clusters_per_resolution <- seurat_object.reactions$clusters_per_resolution[r]

    for(nsid in module_environments$cluster_resolution_pickers$ns)
      updateSelectInput(session=session, inputId=nsid, selected=input$cluster_resolution_picker)
  })
}
