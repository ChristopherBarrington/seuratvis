
10^(0:9) -> major_breaks_log10
(2:9) * rep(major_breaks_log10, each=8) -> minor_breaks_log10

module_environments <- new.env()
module_servers_to_call <- new.env()
seurat_object.reactions <- reactiveValues()
filtered_cells.reactions <- reactiveValues()
seurat_configuration.reactions <- reactiveValues()

#' Reactive lists for dataset filetering
#' 
filtering_parameters.reactions <- reactiveValues()
filtering_arguments.reactions <- reactiveValues()

#' Reactive list of summary statistics of loaded seurat object
#! TODO initialise this when an object is loaded
reference_metrics.rv <- reactiveValues()

#' Plotting reactive values
#' 
plotting_options.rv <- reactiveValues()

#' UI elements by module
#' 
ui_element_ids.env <- new.env()

#' Reactive list of UI element selections
#' 
#' @details
#' Elements of the list can be accessed with \code{session$ns('feature')}.
#' 
#' @param `-picked_feature` name of picked feature from a tab
#' 
selections.rv <- reactiveValues()

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
  session$onSessionEnded(function() {
    message('### session ended')})
}
