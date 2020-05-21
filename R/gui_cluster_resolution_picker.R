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
#'   callModule(cluster_resolution_picker.server, id='page_name')
#' }}
#' 
#' 
#' @rdname cluster_resolution_picker
#' 
cluster_resolution_picker.ui <- function(id, label='Cluster resolutions', include_label_switch=FALSE) {
  message('### cluster_resolution_picker.ui')

  module <- 'cluster_resolution_picker'

  # make unique id for this object
  ns <- NS(namespace=id)
  module_ns <- ns(id=module)

  # create an environment in the seuratvis namespace
  e <- new.env()
  e$id <- id
  assign(x=module_ns, val=e, envir=module_environments)

  module_environments$cluster_resolution_pickers$ns %<>% c(module_ns)
  module_environments$cluster_resolution_pickers$id %<>% c(id)

  # if a label switch is required, make one
  label_switch <- NULL
  if(include_label_switch)
      materialSwitch(inputId=ns(id='label_clusters'), label='Cluster labels', value=TRUE, right=TRUE, status='success') -> label_switch

  # make a drop down selector element
  selectInput(inputId=ns(id='cluster_resolution_picker'), label=label,
              choices='seurat_clusters', selected='seurat_clusters',
              multiple=FALSE) -> dropdown

  # return ui element(s)
  tagList(dropdown, label_switch)
}

#' React to a cluster resolution choice
#' 
#' @details
#' Updates the Seurat object reactive values. Could use \code{SetIdent} instead?
#' 
#' @rdname cluster_resolution_picker
#' 
cluster_resolution_picker.server <- function(input, output, session) {
  message('### cluster_resolution_picker.server')

  # get environments containing variables to run/configure this object
  collect_environments(id=parent.frame()$id, module='cluster_resolution_picker') # provides `seuratvis_env`, `server_env` and `module_env`
  session_server <- get(x='session', env=server_env)

  # react to the cluster resolution selection
  observeEvent(eventExpr=input$cluster_resolution_picker, handlerExpr={
    # create variables for shorthand
    r <- input$cluster_resolution_picker

    # save cluster information in the reactive
    seurat_object.reactions$selected_cluster_resolution <- r
    seurat_object.reactions$selected_clusters_per_resolution <- seurat_object.reactions$clusters_per_resolution[r]

    # update other cluster resolution pickers
    for(nsid in module_environments$cluster_resolution_pickers$ns)
      updateSelectInput(session=session_server, inputId=nsid, selected=input$cluster_resolution_picker)
  })

  # react to the cluster labels switch
  observeEvent(eventExpr=input$label_clusters, handlerExpr={
    seurat_object.reactions$label_clusters <- input$label_clusters
  })

  # update UI when Seurat object is loaded
  observe(x={
    # create variables for shorthand
    seurat <- seurat_object.reactions$seurat
    if(is.null(seurat))
      return(NULL)

    # update the ui element(s)
    cluster_options <- c('seurat_clusters', str_subset(colnames(seurat@meta.data), '_snn_res.'))
    updateSelectInput(session=session, inputId='cluster_resolution_picker',
                      choices=cluster_options, selected='seurat_clusters')})
}
