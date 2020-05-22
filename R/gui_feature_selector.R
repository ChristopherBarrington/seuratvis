#' Select a feature from a dataset
#' 
#' Provides a method to select genes/features of a dataset
#' 
#' @param id unique name of the element
#' @param label text label of the element
#' 
#' @examples
#' 
#' \dontrun{
#' ui <- fluidPage(feature_picker.ui(id='page_name'))
#' server <- function(input, output, session) {
#'   callModule(feature_picker.server, id='page_name')
#' }}
#' 
#' 
#' @rdname feature_picker
#' 
feature_picker.ui <- function(id, label='Feature selection', include_metadata_switch=TRUE) {
  message('### feature_picker.ui')

  module <- 'feature_picker'

  # make unique id for this object
  ns <- NS(namespace=id)
  module_ns <- ns(id=module)

  # create an environment in the seuratvis namespace
  e <- new.env()
  e$id <- id
  assign(x=module_ns, val=e, envir=module_environments)

  module_environments$feature_pickers$ns %<>% c(module_ns)
  module_environments$feature_pickers$id %<>% c(id)

  # make ui elements
  ## if a label switch is required, make one
  metadata_switch <- NULL
  if(include_metadata_switch)
      materialSwitch(inputId=ns(id='list_metadata'), label='Show metadata', value=FALSE, right=FALSE, status='success') -> metadata_switch

  ## feature names autocomplete box
  autocomplete_input(id=ns(id='feature_picker'), label=label,
                     options=NULL, value=NULL,
                     placeholder='Feature') -> dropdown

  ## slider to limit colour range
  sliderInput(inputId=ns(id='value_range'), label='Colour range limits',
              min=0, max=1, step=0.1, value=c(-Inf,Inf)) -> value_range

  # return ui element(s)
  tagList(dropdown, metadata_switch, value_range)
}

#' React to a feature choice
#' 
#' @details
#' Updates the Seurat object reactive values.
#' 
#' @rdname feature_picker
#' 
feature_picker.server <- function(input, output, session) {
  message('### feature_picker.server')

  # get environments containing variables to run/configure this object
  collect_environments(id=parent.frame()$id, module='feature_picker') # provides `seuratvis_env`, `server_env` and `module_env`
  session_server <- get(x='session', env=server_env)

  # react to the feature selection
  observeEvent(eventExpr=input$feature_picker, handlerExpr={
    message('### feature_picker.server-observeEvent-input$feature_picker')

    # create variables for shorthand
    picked <- input$feature_picker
    seurat <- seurat_object.reactions$seurat

    if(is.null(seurat))
      return(NULL)

    # get the values for the selected feature from the loaded Seurat
    picked_feature_values <- FetchData(object=seurat, vars=picked) %>% set_names('value')

    # update the ui element(s)
    ## slider to limit colour range
    max_value <- 1
    if(class(picked_feature_values$value)=='numeric')
      max_value <- max(picked_feature_values$value) %>% add(0.05) %>% round(digits=1)
    updateSliderInput(session=session, inputId='value_range',
                      max=max_value, value=c(-Inf,Inf))

    # save feature information in the reactive
    seurat_object.reactions$picked_feature <- picked
      seurat_object.reactions$picked_feature_values <- picked_feature_values

    # update other feature pickers
    for(nsid in module_environments$feature_pickers$ns)
      updateSelectInput(session=session_server, inputId=nsid, selected=picked)})

  # react to the colour scale limits
  observeEvent(eventExpr=input$value_range, handlerExpr={
    message('### feature_picker.server-observeEvent-input$value_range')

    # update the reactive
    seurat_object.reactions$value_range_limits <- input$value_range})

  # react to the show metadata switch
  observeEvent(eventExpr=input$list_metadata, handlerExpr={
    message('### feature_picker.server-observeEvent-input$list_metadata')

    # create variables for shorthand
    seurat <- seurat_object.reactions$seurat
    if(is.null(seurat))
      return(NULL)
    
    # get the possible features
    selection <- rownames(seurat) %>% sort()
    if(input$list_metadata)
      selection <- colnames(seurat@meta.data)
    default_value <- sample(x=selection, size=1)
    picked_feature_values <- FetchData(object=seurat, vars=default_value) %>% set_names('value')

    # update the ui element(s)
    ## feature names autocomplete box
    update_autocomplete_input(session=session, id='feature_picker',
                              options=selection, value=default_value)

    ## slider to limit colour range
    max_value <- 1
    if(class(picked_feature_values$value)=='numeric')
      max_value <- max(picked_feature_values$value) %>% add(0.05) %>% round(digits=1)
    updateSliderInput(session=session, inputId='value_range',
                      max=max_value, value=c(-Inf,Inf))

    # update the reactive
    seurat_object.reactions$picked_feature <- default_value
    seurat_object.reactions$picked_feature_values <- picked_feature_values})

  # update UI when Seurat object is loaded
  observeEvent(eventExpr=seurat_object.reactions$seurat, handlerExpr={
    message('### feature_picker.server-observeEvent-seurat_object.reactions$seurat')

    # create variables for shorthand
    seurat <- seurat_object.reactions$seurat

    # get the possible features
    selection <- rownames(seurat) %>% sort()
    if(input$list_metadata)
      selection <- colnames(seurat@meta.data) %>% str_subset
    default_value <- sample(x=selection, size=1)
    picked_feature_values <- FetchData(object=seurat, vars=default_value) %>% set_names('value')

    # update the ui element(s)
    ## feature names autocomplete box
    update_autocomplete_input(session=session, id='feature_picker',
                              options=selection, value=default_value)
    
    ## slider to limit colour range
    max_value <- 1
    if(class(picked_feature_values$value)=='numeric')
      max_value <- max(picked_feature_values$value) %>% add(0.05) %>% round(digits=1)
    updateSliderInput(session=session, inputId='value_range',
                      max=max_value, value=c(-Inf,Inf))

    # update the reactive
    seurat_object.reactions$picked_feature <- default_value
    seurat_object.reactions$picked_feature_values <- picked_feature_values})
}
