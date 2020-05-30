
10^(0:9) -> major_breaks_log10
(2:9) * rep(major_breaks_log10, each=8) -> minor_breaks_log10

module_environments <- new.env()
module_servers_to_call <- new.env()
seurat_object.reactions <- reactiveValues()
filtering_parameters.reactions <- reactiveValues()
filtered_cells.reactions <- reactiveValues()
seurat_configuration.reactions <- reactiveValues()

shinyAppServer <- function(input, output, session) {

  # ###############################################################################################
  # scour the session for Seurat objects and populate the UI --------------------------------------
  available_seurat_objects <- find_seurat_objects() # has to be here to access the global environment ... ?

  # ###############################################################################################
  # reactions to tab selection
  ## react to opening tab with a filtered object loaded
  observeEvent(input$sidebarmenu, {
    if(!is.null(seurat_object.reactions$seurat) & input$sidebarmenu=='cell_filtering-tab' && (!is.null(seurat_object.reactions$seurat@misc$cells_filtered) && seurat_object.reactions$seurat@misc$cells_filtered))
      sendSweetAlert(session=session, type='success', html=TRUE,
                     title='Notice', btn_labels='Great!',
                     text=tags$span('It looks like low-quality cells have already been removed from this Seurat object:', tags$h5(tags$code('@misc$cells_filtered == TRUE'))),
                     closeOnClickOutside=TRUE, showCloseButton=FALSE)})

  # ###############################################################################################
  # call servers for modules
  ## modules listed in the module_servers_to_call environment
  module_servers_to_call %<>% as.list()
  for(id in names(module_servers_to_call))
    for(server in module_servers_to_call[[id]])
      callModule(module=get(x=server), id=id)

  ## load the filter_seurat module
  callModule(module=cell_filtering.server, id='seuratvis')

  # ###############################################################################################
  # any code to exectue when the session ends
  session$onSessionEnded(function() {message('### session ended')})
}
