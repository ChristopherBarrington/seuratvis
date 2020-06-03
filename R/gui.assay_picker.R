#' Select assay
#' 
#' Provides a way to select an assay to use
#' 
#' @param id unique name of the element
#' @param label text label of the element
#' 
#' @examples
#' 
#' \dontrun{
#' ui <- fluidPage(assay_picker.ui(id='page_name'))
#' server <- function(input, output, session) {
#'   callModule(assay_picker.server, id='page_name')
#' }}
#' 
#' @rdname assay_picker
#' 
assay_picker.ui <- function(id, label='Assay') {
  sprintf(fmt='### %s-assay_picker.ui', id) %>% message()

  module <- 'assay_picker'

  # make unique id for this object
  ns <- NS(namespace=id)
  module_ns <- ns(id=module)

  # create an environment in the seuratvis namespace
  e <- new.env()
  e$id <- id
  assign(x=module_ns, val=e, envir=module_environments)

  # record the server(s) to call
  module_environments$assay_pickers$ns %<>% c(module_ns) # keep track so other modules can update
  get0(env=module_servers_to_call, x=id) %>% append(sprintf(fmt='%s.server', module)) %>% assign(env=module_servers_to_call, x=id)

  # make a drop down selector element
  selectInput(inputId=ns(id='assay_picker'), label=label,
              choices=NULL, selected=NULL,
              multiple=FALSE) -> dropdown
  
  # return ui element(s)
  dropdown
}

#' React to an assay choice
#' 
#' @details
#' Updates the Seurat object reactive values.
#' 
#' @rdname assay_picker
#' 
assay_picker.server <- function(input, output, session) {
  session$ns('') %>% sprintf(fmt='### %sassay_picker.server') %>% message()

  # get environments containing variables to run/configure this object
  collect_environments(id=parent.frame()$id, module='assay_picker') # provides `seuratvis_env`, `server_env` and `module_env`
  session_server <- get(x='session', env=server_env)

  # react to the reduction method selection
  observeEvent(eventExpr=input$assay_picker, handlerExpr={
    # send a message
    session$ns('') %>% sprintf(fmt='### %sassay_picker.server-observeEvent-input$assay_picker [%s]', input$assay_picker) %>% message('')

    # create variables for shorthand
    assay <- input$assay_picker
    seurat_object.reactions$selected_assay <- assay

    # update other reduction method pickers
    for(nsid in module_environments$assay_pickers$ns)
      updateSelectInput(session=session_server, inputId=nsid, selected=assay)

    # set the default assay for the reactive seurat object
    if(assay!='' && !is.null(seurat_object.reactions$seurat))
      DefaultAssay(seurat_object.reactions$seurat) <- assay

    # update the reactive
    seurat_object.reactions$reference_metrics$n_features <- nrow(seurat_object.reactions$seurat)})

  # update UI when Seurat object is loaded
  observeEvent(eventExpr=seurat_object.reactions$seurat, handlerExpr={
    # send a message
    session$ns('') %>% sprintf(fmt='### %sassay_picker.server-observeEvent-seurat_object.reactions$seurat [%s]', seurat_object.reactions$formatted.project.name) %>% message()

    # create variables for shorthand
    seurat <- seurat_object.reactions$seurat

    # update the ui element(s)
    updateSelectInput(session=session, inputId='assay_picker',
                      choices=Assays(seurat), selected=DefaultAssay(seurat))})
}
