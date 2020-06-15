
10^(0:9) -> major_breaks_log10
(2:9) * rep(major_breaks_log10, each=8) -> minor_breaks_log10

module_environments <- new.env()
module_servers_to_call <- new.env()
seurat_object.reactions <- reactiveValues()

#' UI elements by module
#' 
ui_element_ids.env <- new.env()

shinyAppServer <- function(input, output, session) {

  # ###############################################################################################
  # scour the session for Seurat objects and populate the UI --------------------------------------
  #! TODO: this is not the right place for this. maybe a call to utils::globalVariables() in seuratvis() ... ?
  available_seurat_objects <- find_seurat_objects() # has to be here to access the global environment ... ?

  # ###############################################################################################
  # call servers for modules
  ## special case to load the seurat handler tab's server first
  callModule(module=available_seurats.server, id='load_dataset')
  seurat <- callModule(module=seurat_object.server, id='load_dataset')  # callModule(module=load_a_seurat.server, id='load_dataset')

  ## load the filter_seurat module
  cell_filtering <- callModule(module=cell_filtering.server, id='seuratvis', seurat=seurat)

  ## modules listed in the module_servers_to_call environment
  module_servers_to_call %<>% as.list()
  for(id in names(module_servers_to_call))
    for(server in module_servers_to_call[[id]]) {
      if(TRUE) print(server)
      callModule(module=get(x=server), id=id, seurat=seurat, cell_filtering=cell_filtering)
    }

  # ###############################################################################################
  # reactions to tab selection
  ## react to opening tab with a filtered object loaded
  observeEvent(input$sidebarmenu, {
    if(!is.null(seurat$object) & input$sidebarmenu=='cell_filtering-tab' && (!is.null(seurat$object@misc$cells_filtered) && seurat$object@misc$cells_filtered))
      sendSweetAlert(session=session, type='success', html=TRUE,
                     title='Notice', btn_labels='Great!',
                     text=tags$span('It looks like low-quality cells have already been removed from this Seurat object:', tags$h5(tags$code('@misc$cells_filtered == TRUE'))),
                     closeOnClickOutside=TRUE, showCloseButton=FALSE)})

  # ###############################################################################################
  # any code to exectue when the session ends
  session$onSessionEnded(function() {
    message('### session ended')})
}
