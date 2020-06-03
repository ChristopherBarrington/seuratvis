#' Select reduction method
#' 
#' Provides a way to select a reduction method to plot
#' 
#' @param id unique name of the element
#' @param label text label of the element
#' 
#' @examples
#' 
#' \dontrun{
#' ui <- fluidPage(reduction_method_picker.ui(id='page_name'))
#' server <- function(input, output, session) {
#'   callModule(reduction_method_picker.server, id='page_name')
#' }}
#'
#' @rdname reduction_method_picker
#' 
reduction_method_picker.ui <- function(id, label='Reduction method') {
  sprintf(fmt='### %s-reduction_method_picker.ui', id) %>% message()

  module <- 'reduction_method_picker'

  # make unique id for this object
  ns <- NS(namespace=id)
  module_ns <- ns(id=module)

  # create an environment in the seuratvis namespace
  e <- new.env()
  e$id <- id
  assign(x=module_ns, val=e, envir=module_environments)

  # record the server(s) to call
  module_environments$reduction_method_pickers$ns %<>% c(module_ns) # keep track so other modules can update
  get0(env=module_servers_to_call, x=id) %>% append(sprintf(fmt='%s.server', module)) %>% assign(env=module_servers_to_call, x=id)

  # make a drop down selector element
  selectInput(inputId=ns(id='reduction_method_picker'), label=label,
              choices=NULL, selected=NULL,
              multiple=FALSE) -> dropdown
  
  # return ui element(s)
  dropdown
}

#' React to a reduction method choice
#' 
#' @details
#' Updates the Seurat object reactive values.
#' 
#' @rdname reduction_method_picker
#' 
reduction_method_picker.server <- function(input, output, session) {
  session$ns('') %>% sprintf(fmt='### %sreduction_method_picker.server') %>% message()

  # get environments containing variables to run/configure this object
  collect_environments(id=parent.frame()$id, module='reduction_method_picker') # provides `seuratvis_env`, `server_env` and `module_env`
  session_server <- get(x='session', env=server_env)

  # react to the reduction method selection
  observeEvent(eventExpr=input$reduction_method_picker, handlerExpr={
    message('### reduction_method_picker.server-observeEvent-input$reduction_method_picker')

    # create varaibles for shorthand
    dimred_method <- input$reduction_method_picker
    seurat_object.reactions$selected_reduction_method <- dimred_method

    # update other reduction method pickers
    for(nsid in module_environments$reduction_method_pickers$ns)
      updateSelectInput(session=session_server, inputId=nsid, selected=dimred_method)

    # pull out the reduction and put a data.frame in the seurat object reactions
    seurat <- seurat_object.reactions$seurat
    if(dimred_method!='' && !is.null(seurat))
      seurat@reductions[[dimred_method]]@cell.embeddings[,1:2] %>%
        as.data.frame() %>%
        set_names(c('DIMRED_1','DIMRED_2')) %>%
        cbind(seurat@meta.data) -> seurat_object.reactions$dimred})

  # update UI when Seurat object is loaded
  observeEvent(eventExpr=seurat_object.reactions$seurat, handlerExpr={
    sprintf(fmt='### reduction_method_picker.server-observeEvent-seurat_object.reactions$seurat [%s]', seurat_object.reactions$formatted.project.name) %>% message()

    # create varaibles for shorthand
    seurat <- seurat_object.reactions$seurat
    reductions <- Reductions(seurat)

    if(length(reductions)==0)
      return(NULL)

    # update the ui element(s)
    updateSelectInput(session=session, inputId='reduction_method_picker',
                      choices=Reductions(seurat),
                      selected=Seurat:::DefaultDimReduc(seurat))

    # initialise the dimension reduced map
    seurat@reductions[[Seurat:::DefaultDimReduc(seurat)]]@cell.embeddings[,1:2] %>%
        as.data.frame() %>%
        set_names(c('DIMRED_1','DIMRED_2')) %>%
        cbind(seurat@meta.data) -> seurat_object.reactions$dimred
  })
}
