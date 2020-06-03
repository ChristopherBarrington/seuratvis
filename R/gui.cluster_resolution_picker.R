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
#' @rdname cluster_resolution_picker
#' 
cluster_resolution_picker.ui <- function(id, label='Cluster resolutions', include_label_switch=FALSE) {
  sprintf(fmt='### %s-cluster_resolution_picker.ui', id) %>% message()

  module <- 'cluster_resolution_picker'

  # make unique id for this object
  ns <- NS(namespace=id)
  module_ns <- ns(id=module)

  # create an environment in the seuratvis namespace
  e <- new.env()
  e$id <- id
  assign(x=module_ns, val=e, envir=module_environments)

  # record the server(s) to call
  module_environments$cluster_resolution_pickers$ns %<>% c(module_ns) # keep track so other modules can update
  get0(env=module_servers_to_call, x=id) %>% append(sprintf(fmt='%s.server', module)) %>% assign(env=module_servers_to_call, x=id)

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
  session$ns('') %>% sprintf(fmt='### %scluster_resolution_picker.server') %>% message()

  # get environments containing variables to run/configure this object
  collect_environments(id=parent.frame()$id, module='cluster_resolution_picker') # provides `seuratvis_env`, `server_env` and `module_env`
  session_server <- get(x='session', env=server_env)

  # react to the cluster resolution selection
  observeEvent(eventExpr=input$cluster_resolution_picker, handlerExpr={
    # make sure these elements are defined
    req(seurat_object.reactions$seurat)

    # send a message
    session$ns('') %>% sprintf(fmt='### %scluster_resolution_picker.server-observeEvent-input$cluster_resolution_picker [%s]', input$cluster_resolution_picker) %>% message('')

    # create variables for shorthand
    r <- input$cluster_resolution_picker
    seurat <- seurat_object.reactions$seurat

    # save cluster information in the reactive
    seurat_object.reactions$selected_cluster_resolution <- r
    seurat_object.reactions$selected_clusters_per_resolution <- seurat_object.reactions$clusters_per_resolution[r]
    seurat_object.reactions$picked_cluster_resolution_idents <- FetchData(object=seurat, vars=r) %>% set_names('ident')

    # update other cluster resolution pickers
    for(nsid in module_environments$cluster_resolution_pickers$ns)
      updateSelectInput(session=session_server, inputId=nsid, selected=input$cluster_resolution_picker)})

  # react to the cluster labels switch
  observeEvent(eventExpr=input$label_clusters, handlerExpr={
    # send a message
    session$ns('') %>% sprintf(fmt='### %scluster_resolution_picker.server-observeEvent-input$label_clusters [%s]', input$label_clusters) %>% message('')

    seurat_object.reactions$label_clusters <- input$label_clusters})

  # update UI when Seurat object is loaded
  observeEvent(eventExpr=seurat_object.reactions$seurat, handlerExpr={
    # send a message
    session$ns('') %>% sprintf(fmt='### %sreduction_method_picker.server-observeEvent-seurat_object.reactions$seurat [%s]', seurat_object.reactions$formatted.project.name) %>% message()

    # create variables for shorthand
    seurat <- seurat_object.reactions$seurat

    # update the ui element(s)
    cluster_options <- c('seurat_clusters', str_subset(colnames(seurat@meta.data), '_snn_res.'))
    updateSelectInput(session=session, inputId='cluster_resolution_picker',
                      choices=cluster_options, selected='seurat_clusters')

    # count the clusters per resolution and initialise reactive
    select_at(seurat@meta.data, vars(all_of(cluster_options))) %>%
      mutate_all(function(x) {as.character(x) %>% as.numeric()}) %>%
      gather(key='cluster_set', value='ID') %>%
      group_by(cluster_set) %>%
      summarise(N=length(unique(ID))) %>%
      deframe() -> clusters_per_resolution

    seurat_object.reactions$clusters_per_resolution <- clusters_per_resolution
    seurat_object.reactions$selected_clusters_per_resolution <- clusters_per_resolution['seurat_clusters']
    seurat_object.reactions$picked_cluster_resolution_idents <- FetchData(object=seurat, vars='seurat_clusters') %>% set_names('ident')})
}
