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
  sprintf(fmt='### %s-feature_picker.ui', id) %>% message()

  module <- 'feature_picker'

  # make unique id for this object
  ns <- NS(namespace=id)
  module_ns <- ns(id=module)

  # create an environment in the seuratvis namespace
  e <- new.env()
  e$id <- id
  e$include_metadata_switch <- include_metadata_switch
  assign(x=module_ns, val=e, envir=module_environments)

  # track the re-used UI elements in each namespace
  get0(env=ui_element_ids.env, x=NS(namespace=module, id='feature_picker_feature_names')) %>%
    append(ns(id='feature_picker_feature_names')) %>%
    assign(env=ui_element_ids.env, x=NS(namespace=module, id='feature_picker_feature_names'))

  # record the server(s) to call
  get0(env=module_servers_to_call, x=id) %>% append(sprintf(fmt='%s.server', module)) %>% assign(env=module_servers_to_call, x=id)

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
  session$ns('') %>% sprintf(fmt='### %sfeature_picker.server') %>% message()

  # get environments containing variables to run/configure this object
  collect_environments(id=parent.frame()$id, module='feature_picker') # provides `seuratvis_env`, `server_env` and `module_env`
  session_server <- get(x='session', env=server_env)

  # react to the feature selection
  ## if a feature is selected, copy it to the reactive
  observeEvent(eventExpr={input$feature_picker_feature_names }, handlerExpr={
    # make sure these elements are defined
    req(input$feature_picker_feature_names)
    
    # send a message
    sprintf(fmt='### %sfeature_picker.server-observeEvent-input$feature_picker_feature_names [%s]', session$ns(''), input$feature_picker_feature_names) %>% message()
    
    if(is.null(input$list_metadata) || !input$list_metadata)
      selections.rv[[session$ns('picked_feature')]] <- input$feature_picker_feature_names})

  ## if a metadata column is selected, copy it to the reactive
  observeEvent(eventExpr=input$feature_picker_metadata, handlerExpr={
    # make sure these elements are defined
    req(input$feature_picker_metadata)

    # send a message
    sprintf(fmt='### %sfeature_picker.server-observeEvent-input$feature_picker_metadata [%s]', session$ns(''), input$feature_picker_metadata) %>% message()
    
    if(!is.null(input$list_metadata) && input$list_metadata)
      selections.rv[[session$ns('picked_feature')]] <- input$feature_picker_metadata})

  ## if the metadata switch is toggled, set the picked feature
  observeEvent(eventExpr=input$list_metadata, handlerExpr={
    # make sure these elements are defined
    req(input$list_metadata)
    req(selections.rv[[session$ns('picked_feature_previous')]])

    # send a message
    sprintf(fmt='### %sfeature_picker.server-observeEvent-input$list_metadata [%s]', session$ns(''), input$list_metadata) %>% message()

    # pick the feature to revert to
    ## if metadata switch is true, get the value of the metadata dropdown
    ## if metadata switch is false, get the previously shown feature (autocomplete_input return empty in this case)
    picked_feature <- ifelse(input$list_metadata, input$feature_picker_metadata, selections.rv[[session$ns('picked_feature_previous')]])
    selections.rv[[session$ns('picked_feature_previous')]] <- selections.rv[[session$ns('picked_feature_previous')]]

    # update the reactive
    selections.rv[[session$ns('picked_feature')]] <- picked_feature})

  ## use the selected feature (it may be a feature or metadata)
  observeEvent(eventExpr=selections.rv[[session$ns('picked_feature')]], handlerExpr={
    # make sure these elements are defined
    req(seurat_object.reactions$seurat)
    req(selections.rv[[session$ns('picked_feature')]])

    # send a message
    sprintf(fmt='### %sfeature_picker.server-observeEvent-seurat_object.reactions$picked_feature [%s]', session$ns(''), selections.rv[[session$ns('picked_feature_previous')]]) %>% message()

    # create variables for shorthand
    picked <- selections.rv[[session$ns('picked_feature')]]
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
    selections.rv[[session$ns('picked_feature_values')]] <- picked_feature_values
    selections.rv[[session$ns('value_range_limits')]] <- c(min_value, max_value)})

  # react to the colour scale limits
  observeEvent(eventExpr=input$value_range, handlerExpr={
    # make sure these elements are defined
    req(selections.rv[[session$ns('value_range_limits')]])

    # send a message
    sprintf(fmt='### %sfeature_picker.server-observeEvent-input$value_range [%s]', session$ns(''), str_c(input$value_range, collapse=',')) %>% message()

    # update the reactive
    if(!identical(seurat_object.reactions$value_range_limits, input$value_range))
      selections.rv[[session$ns('value_range_limits')]] <- input$value_range})

  # update UI when Seurat object is loaded
  observeEvent(eventExpr=seurat_object.reactions$seurat, handlerExpr={
    # send a message
    sprintf(fmt='### %sfeature_picker.server-observeEvent-seurat_object.reactions$seurat [%s]', session$ns(''), seurat_object.reactions$formatted.project.name) %>% message()

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
    seurat_object.reactions$feature_picker_features <- feature_picker_options$features
    seurat_object.reactions$feature_picker_metadata <- feature_picker_options$metadata

    selections.rv[[session$ns('picked_feature')]] <- picked_feature
    selections.rv[[session$ns('picked_feature_previous')]] <- picked_feature})
}
