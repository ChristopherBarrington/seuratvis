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
cluster_resolution_picker.server <- function(input, output, session, seurat, ...) {
  session$ns('') %>% sprintf(fmt='### %scluster_resolution_picker.server') %>% message()

  # get environments containing variables to run/configure this object
  collect_environments(id=parent.frame()$id, module='cluster_resolution_picker') # provides `seuratvis_env`, `server_env` and `module_env`
  session_server <- get(x='session', env=server_env)

  # react to the cluster resolution selection
  observeEvent(eventExpr=input$cluster_resolution_picker, handlerExpr={
    # make sure these elements are defined
    req(seurat$object)

    # send a message
    session$ns('') %>% sprintf(fmt='### %scluster_resolution_picker.server-observeEvent-input$cluster_resolution_picker [%s]', input$cluster_resolution_picker) %>% message()

    # create variables for shorthand
    r <- input$cluster_resolution_picker
    object <- seurat$object

    # save cluster information in the reactive
    seurat$selected_clusters_per_resolution <- seurat$clusters_per_resolution[r]
    seurat$picked_cluster_resolution_idents <- FetchData(object=object, vars=r) %>% set_names('ident')})

  # update UI when Seurat object is loaded
  observeEvent(eventExpr=seurat$object, handlerExpr={
    # send a message
    session$ns('') %>% sprintf(fmt='### %scluster_resolution_picker.server-observeEvent-seurat$object [%s]', seurat$formatted_project) %>% message()

    # create variables for shorthand
    object <- seurat$object

    # update the ui element(s)
    cluster_options <- c('seurat_clusters', str_subset(colnames(object@meta.data), '_snn_res.'))
    updateSelectInput(session=session, inputId='cluster_resolution_picker', choices='-spoof-')
    updateSelectInput(session=session, inputId='cluster_resolution_picker',
                      choices=cluster_options, selected='seurat_clusters')

    # count the clusters per resolution and initialise reactive
    select_at(object@meta.data, vars(all_of(cluster_options))) %>%
      mutate_all(function(x) {as.character(x) %>% as.numeric()}) %>%
      gather(key='cluster_set', value='ID') %>%
      group_by(cluster_set) %>%
      summarise(N=length(unique(ID))) %>%
      deframe() -> clusters_per_resolution

    seurat$clusters_per_resolution <- clusters_per_resolution
    seurat$selected_clusters_per_resolution <- clusters_per_resolution['seurat_clusters']
    seurat$picked_cluster_resolution_idents <- FetchData(object=object, vars='seurat_clusters') %>% set_names('ident')})
}
