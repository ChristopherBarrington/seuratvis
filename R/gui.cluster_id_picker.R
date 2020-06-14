#' Pick clusters from selected resolution
#' 
#' Provides a method to select clusters(s)
#' 
#' @param id unique name of the element
#' @param label text label of the element
#' 
#' @examples
#' 
#' \dontrun{
#' ui <- fluidPage(cluster_id_picker.ui(id='page_name'))
#' server <- function(input, output, session) {
#'   callModule(cluster_id_picker.server, id='page_name')
#' }}
#' 
#' @rdname cluster_id_picker
#' 
cluster_id_picker.ui <- function(id, label='Cluster selection') {
  sprintf(fmt='### %s-cluster_id_picker.ui', id) %>% message()

  module <- 'cluster_id_picker'

  # make unique id for this object
  ns <- NS(namespace=id)
  module_ns <- ns(id=module)

  # create an environment in the seuratvis namespace
  e <- new.env()
  e$id <- id
  assign(x=module_ns, val=e, envir=module_environments)

  # record the server(s) to call
  get0(env=module_servers_to_call, x=id) %>% append(sprintf(fmt='%s.server', module)) %>% assign(env=module_servers_to_call, x=id)

  # make ui elements
  pickerInput(inputId=ns(id='cluster_id_picker'), label=label, choices=NULL, multiple=TRUE,
              options=list(`actions-box`=TRUE, header='Cluster selection', title='Cluster selection',
                           `selected-text-format`='count>5', `count-selected-text`='{0} cluster(s)')) -> input

  # return ui element(s)
  tagList(input)
}

#' React to cluster selection
#' 
#' @details
#' Updates reactive values.
#' 
#' @import gtools
#' 
#' @rdname cluster_id_picker
#' 
cluster_id_picker.server <- function(input, output, session) {
  session$ns('') %>% sprintf(fmt='### %scluster_id_picker.server') %>% message()

  # get environments containing variables to run/configure this object
  collect_environments(id=parent.frame()$id, module='cluster_id_picker') # provides `seuratvis_env`, `server_env` and `module_env`
  session_server <- get(x='session', env=server_env)

  # update UI when Seurat object is loaded
  observeEvent(eventExpr=seurat_object.reactions$picked_cluster_resolution_idents, handlerExpr={
    # send a message
    session$ns('') %>% sprintf(fmt='### %scluster_id_picker.server-observeEvent-seurat_object.reactions$picked_cluster_resolution_idents [%s]', seurat_object.reactions$formatted.project.name) %>% message('')

    # create variables for shorthand
    idents <- seurat_object.reactions$picked_cluster_resolution_idents %>% pluck('ident') %>% levels() %>% mixedsort()

    # update the ui
    updatePickerInput(session=session, inputId='cluster_id_picker', choices=idents, selected=idents)})
}
