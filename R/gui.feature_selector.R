#' Select a feature from a dataset
#' 
#' Provides a method to select genes/features/metadata of a dataset
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
  e$include_metadata_switch <- include_metadata_switch
  assign(x=module_ns, val=e, envir=module_environments)

  module_environments$feature_pickers$ns %<>% c(module_ns)
  module_environments$feature_pickers$id %<>% c(id)

  # make ui elements
  ## if a label switch is required, make one
  metadata_switch <- NULL
  if(include_metadata_switch)
      materialSwitch(inputId=ns(id='list_metadata'), label='Show metadata', value=FALSE, right=FALSE, status='success') -> metadata_switch

  ## feature names autocomplete box
  autocomplete_input(id=ns(id='feature_picker_feature_names'), label='Feature name', options=NULL, value=NULL, placeholder='Feature') %>%
    conditionalPanel(condition=sprintf('!input["%s"]', ns(id='list_metadata'))) -> feature_names_picker

  ## metadata names drop down box
  selectizeInput(inputId=ns(id='feature_picker_metadata'), label='Metadata', choices=NULL, selected=NULL) %>%
    conditionalPanel(condition=sprintf('input["%s"]', ns(id='list_metadata'))) -> metadata_picker

  ## slider to limit colour range
  sliderInput(inputId=ns(id='value_range'), label='Colour range limits',
              min=0, max=1, step=0.1, value=c(-Inf,Inf)) -> value_range

  # return ui element(s)
  # tagList(dropdown, metadata_switch, value_range)
  tagList(feature_names_picker,
          metadata_picker,
          metadata_switch,
          value_range)
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
  ## if a feature is selected, copy it to the reactive
  observeEvent(eventExpr=input$feature_picker_feature_names, handlerExpr={
    sprintf(fmt='### feature_picker.server-observeEvent-input$feature_picker_feature_names', input$feature_picker_feature_names) %>% message()
    if(is.null(input$list_metadata) || !input$list_metadata)
      seurat_object.reactions$picked_feature <- input$feature_picker_feature_names
    })

  ## if a metadata column is selected, copy it to the reactive
  observeEvent(eventExpr=input$feature_picker_metadata, handlerExpr={
    sprintf(fmt='### feature_picker.server-observeEvent-input$feature_picker_metadata [%s]', input$feature_picker_metadata) %>% message()
    if(!is.null(input$list_metadata) && input$list_metadata)
      seurat_object.reactions$picked_feature <- input$feature_picker_metadata})

  ## if the metadata switch is toggled, set the picked feature
  observeEvent(eventExpr=input$list_metadata, handlerExpr={
    sprintf(fmt='### feature_picker.server-observeEvent-input$list_metadata [%s]', input$list_metadata) %>% message()

    # make sure these elements are defined
    req(seurat_object.reactions$picked_feature_previous)

    # pick the feature to revert to
    ## if metadata switch is true, get the value of the metadata dropdown
    ## if metadata switch is false, get the previously shown feature (autocomplete_input return empty in this case)
    picked_feature <- ifelse(input$list_metadata, input$feature_picker_metadata, seurat_object.reactions$picked_feature_previous)
    seurat_object.reactions$picked_feature_previous <- seurat_object.reactions$picked_feature

    # update the reactive
    seurat_object.reactions$picked_feature <- picked_feature})

  ## use the selected feature (it may be a feature or metadata)
  observeEvent(eventExpr=seurat_object.reactions$picked_feature, handlerExpr={  
    sprintf(fmt='### feature_picker.server-observeEvent-seurat_object.reactions$picked_feature [%s]', seurat_object.reactions$picked_feature) %>% message()

    # make sure these elements are defined
    req(seurat_object.reactions$seurat)

    # create variables for shorthand
    picked <- seurat_object.reactions$picked_feature
    seurat <- seurat_object.reactions$seurat

    # get the values for the selected feature from the loaded Seurat
    picked_feature_values <- FetchData(object=seurat, vars=picked) %>% set_names('value')

    # update the ui element(s)
    ## slider to limit colour range
    min_value <- 0
    max_value <- 1
    if(class(picked_feature_values$value)=='numeric') {
      min_value <- min(picked_feature_values$value) %>% subtract(0.05) %>% round(digits=1)
      max_value <- max(picked_feature_values$value) %>% add(0.05) %>% round(digits=1)
    }

    updateSliderInput(session=session, inputId='value_range',
                      min=min_value, max=max_value, value=c(-Inf,Inf))

    # save feature information in the reactive
    seurat_object.reactions$picked_feature_values <- picked_feature_values
    seurat_object.reactions$value_range_limits <- c(min_value, max_value)})

  # react to the colour scale limits
  observeEvent(eventExpr=input$value_range, handlerExpr={
    input$value_range %>%
      str_c(collapse=',') %>%
      sprintf(fmt='### feature_picker.server-observeEvent-input$value_range [%s]') %>%
      message()

    # make sure these elements are defined
    req(seurat_object.reactions$value_range_limits)

    # update the reactive
    if(!identical(seurat_object.reactions$value_range_limits, input$value_range))
      seurat_object.reactions$value_range_limits <- input$value_range})

  # update UI when Seurat object is loaded
  observeEvent(eventExpr=seurat_object.reactions$seurat, handlerExpr={
    sprintf(fmt='### feature_picker.server-observeEvent-seurat_object.reactions$seurat [%s]', seurat_object.reactions$formatted.project.name) %>% message()

    # create variables for shorthand
    seurat <- seurat_object.reactions$seurat

    # get the possible features and values
    ## get names of features and metadata
    list(features=rownames(seurat),
         metadata=colnames(seurat@meta.data)) -> feature_picker_options

    ## pick a random feature and metadata column
    feature_picker_options %>% lapply(sample, size=1) -> feature_picker_selected

    ## get the values of the random features
    feature_picker_selected %>% lapply(FetchData, object=seurat) %>% lapply(set_names, nm='value') -> feature_picker_data

    # pick a feature to display: features or metadata
    picked_feature <- feature_picker_selected$features
    picked_feature_values <- feature_picker_data$features
    if(!is.null(input$list_metadata) && input$list_metadata) {
      picked_feature <- feature_picker_selected$metadata
      picked_feature_values <- feature_picker_data$metadata
    }

    # update the ui element(s)
    ## feature names autocomplete box
    update_autocomplete_input(session=session, id='feature_picker_feature_names',
                              options=feature_picker_options$features, value=feature_picker_selected$features)
    
    ## metadata names dropdown box
    updateSelectizeInput(session=session, inputId='feature_picker_metadata',
                         choices=feature_picker_options$metadata, selected=feature_picker_selected$metadata)

    # update the reactive
    seurat_object.reactions$picked_feature <- picked_feature
    seurat_object.reactions$picked_feature_previous <- picked_feature
    seurat_object.reactions$feature_picker_features <- feature_picker_options$features
    seurat_object.reactions$feature_picker_metadata <- feature_picker_options$metadata})
}
