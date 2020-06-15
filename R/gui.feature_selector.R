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
feature_picker.ui <- function(id, label='Feature selection', include_features=TRUE, include_metadata=TRUE, include_values_range=TRUE, features_regex='.*', metadata_regex='.*') {
  sprintf(fmt='### %s-feature_picker.ui', id) %>% message()

  module <- 'feature_picker'

  # make unique id for this object
  ns <- NS(namespace=id)
  module_ns <- ns(id=module)

  # create an environment in the seuratvis namespace
  e <- new.env()
  e$id <- id
  e$include_values_range <- include_values_range
  e$features_regex <- features_regex
  e$metadata_regex <- metadata_regex
  assign(x=module_ns, val=e, envir=module_environments)

  # record the server(s) to call
  get0(env=module_servers_to_call, x=id) %>% append(sprintf(fmt='%s.server', module)) %>% assign(env=module_servers_to_call, x=id)

  # make ui elements
  ## feature names autocomplete box
  autocomplete_input(id=ns(id='feature_picker_feature_names'), label=NULL, options=NULL, value=NULL, placeholder='Feature') %>%
    conditionalPanel(condition=sprintf('input["%s"]=="features"', ns(id='feature_type'))) -> feature_names_picker_conditional

  ## metadata names drop down box
  selectizeInput(inputId=ns(id='feature_picker_metadata'), label=NULL, choices=NULL, selected=NULL, multiple=FALSE) %>%
    conditionalPanel(condition=sprintf('input["%s"]=="metadata"', ns(id='feature_type'))) -> metadata_picker_conditional

  ## gene modules drop down box
  selectizeInput(inputId=ns(id='feature_picker_gene_module'), label=NULL, choices=NULL, selected=NULL, multiple=FALSE) %>%
    conditionalPanel(condition=sprintf('input["%s"]=="gene_modules"', ns(id='feature_type'))) -> gene_module_picker_conditional

  ## slider to limit colour range
  sliderInput(inputId=ns(id='value_range'), label='Colour range limits',
              min=0, max=1, step=0.1, value=c(-Inf,Inf)) -> value_range

  ## checkbox for feature type
  radioGroupButtons(inputId=ns(id='feature_type'), status='primary', width='100%',
                    choices=list(`Features`='features', `Metadata`='metadata', `Gene modules`='gene_modules'),
                    selected='features',
                    checkIcon=list(yes=icon('ok', lib='glyphicon'))) -> feature_type_picker

  ## hidden text box to serve app
  textInput(inputId=ns('picked_feature'), label='picked feature') -> picked_feature_text_input

  # return ui element(s)
  tagList(h6('Feature selection'),
          feature_type_picker,
          feature_names_picker_conditional,
          metadata_picker_conditional,
          gene_module_picker_conditional,
          if(include_values_range) value_range,
          hidden(picked_feature_text_input))
}

#' React to a feature choice
#' 
#' @details
#' Updates the Seurat object reactive values.
#' 
#' @rdname feature_picker
#' 
feature_picker.server <- function(input, output, session, seurat, ...) {
  session$ns('') %>% sprintf(fmt='### %sfeature_picker.server') %>% message()

  # get environments containing variables to run/configure this object
  collect_environments(id=parent.frame()$id, module='feature_picker') # provides `seuratvis_env`, `server_env` and `module_env`
  session_server <- get(x='session', env=server_env)
  previously_picked_feature <- reactiveValues()
  tab <- parent.frame()$id
  seurat$picked_feature_values <- list()

  # react to the feature selection
  ## if a feature is selected, copy it to the reactive
  observeEvent(eventExpr=input$feature_picker_feature_names, handlerExpr={
    # make sure these elements are defined
    req(input$feature_picker_feature_names)
    
    # send a message
    sprintf(fmt='### %sfeature_picker.server-observeEvent-input$feature_picker_feature_names [%s]', session$ns(''), input$feature_picker_feature_names) %>% message()
    
    # update hidden ui element
    if(input$feature_type=='features')
      updateTextInput(session=session, inputId='picked_feature', value=input$feature_picker_feature_names)})

  ## if a metadata column is selected, copy it to the reactive
  observeEvent(eventExpr=input$feature_picker_metadata, handlerExpr={
    # make sure these elements are defined
    req(input$feature_picker_metadata)

    # send a message
    sprintf(fmt='### %sfeature_picker.server-observeEvent-input$feature_picker_metadata [%s]', session$ns(''), input$feature_picker_metadata) %>% message()

    # update hidden ui element
    if(input$feature_type=='metadata')
      updateTextInput(session=session, inputId='picked_feature', value=input$feature_picker_metadata)})

  ## if a gene module column is selected, copy it to the reactive
  observeEvent(eventExpr=input$feature_picker_gene_module, handlerExpr={
    # make sure these elements are defined
    req(input$feature_picker_gene_module)

    # send a message
    sprintf(fmt='### %sfeature_picker.server-observeEvent-input$feature_picker_gene_module [%s]', session$ns(''), input$feature_picker_gene_module) %>% message()

    # update hidden ui element
    if(input$feature_type=='gene_modules')
      updateTextInput(session=session, inputId='picked_feature', value=input$feature_picker_gene_module)})

  ## update the hidden ui element when a feature type is selected
  observeEvent(eventExpr=input$feature_type, ignoreInit=TRUE, handlerExpr={
    # send a message
    sprintf(fmt='### %sfeature_picker.server-observeEvent-input$feature_type [%s]', session$ns(''), input$feature_type) %>% message()

    # pick the feature to revert to
    input_name <- switch(input$feature_type, features='feature_picker_feature_names', metadata='feature_picker_metadata', gene_modules='feature_picker_gene_module')
    picked_feature <- input[[input_name]]
    if(input_name=='feature_picker_feature_names') # autocomplete_input can return empty in this case before it has been typed into
      picked_feature <- previously_picked_feature[[tab]][[input$feature_type]]

    # update hidden ui element
    updateTextInput(session=session, inputId='picked_feature', value=picked_feature)})

  ## use the selected feature (it may be a feature or metadata)
  observeEvent(eventExpr=input$picked_feature, handlerExpr={
    # make sure these elements are defined
    req(seurat$object)
    req(input$feature_type)

    # send a message
    sprintf(fmt='### %sfeature_picker.server-observeEvent-input$picked_feature [%s]', session$ns(''), input$picked_feature) %>% message()

    # create variables for shorthand
    picked <- input$picked_feature

    # track the history
    previously_picked_feature[[tab]][[input$feature_type]] <- picked

    # get the values for the selected feature from the loaded Seurat
    if(input$feature_type=='gene_modules') {
      picked_feature_values <- seurat$gene_modules %>% select(picked) %>% set_names('value')
    } else {
      picked_feature_values <- FetchData(object=seurat$object, vars=picked) %>% set_names('value')
    }

    # save feature information in the reactive
    #! TODO: this invalidates the features per cluster boxplot but not the umap
    seurat$picked_feature_values[[tab]] <- picked_feature_values

    # update the ui element(s)
    ## slider to limit colour range
    #! TODO: this invalidates the umaps but not the features per cluster boxplot
    min_value <- 0
    max_value <- 1
    if(class(picked_feature_values$value)=='numeric') {
      min_value <- min(picked_feature_values$value) %>% subtract(0.05) %>% round(digits=1)
      max_value <- max(picked_feature_values$value) %>% add(0.05) %>% round(digits=1)
    }

    updateSliderInput(session=session, inputId='value_range',
                      min=min_value, max=max_value, value=c(-Inf,Inf))})

  # update UI when Seurat object is loaded
  observeEvent(eventExpr=seurat$object, handlerExpr={
    # send a message
    sprintf(fmt='### %sfeature_picker.server-observeEvent-seurat$object [%s]', session$ns(''), seurat$formatted_project) %>% message()

    # get the possible features and values
    ## get names of features and metadata
    list(features=rownames(seurat$object),
         metadata=colnames(seurat$metadata),
         gene_modules=colnames(seurat$gene_modules)) -> feature_picker_options

    ## filter the options using the regex
    feature_picker_options$features %<>% str_subset(pattern=regex(pattern=module_env$features_regex, ignore_case=TRUE))
    feature_picker_options$metadata %<>% str_subset(pattern=regex(pattern=module_env$metadata_regex, ignore_case=TRUE))

    ## pick a random feature and metadata and gene module column
    feature_picker_options %>%
      lapply(sample, size=1) -> feature_picker_selected

    # pick a feature to display: features or metadata or gene_module
    picked_feature <- feature_picker_selected %>% pluck(input$feature_type)

    # update the ui element(s)
    ## feature names autocomplete box
    update_autocomplete_input(session=session, id='feature_picker_feature_names', options='-spoof-')
    update_autocomplete_input(session=session, id='feature_picker_feature_names',
                              options=feature_picker_options$features, value=feature_picker_selected$features)
    
    ## metadata names dropdown box
    updateSelectizeInput(session=session, inputId='feature_picker_metadata', choices='-spoof-')
    updateSelectizeInput(session=session, inputId='feature_picker_metadata',
                         choices=feature_picker_options$metadata, selected=feature_picker_selected$metadata)

    ## gene modules dropdown box
    updateSelectizeInput(session=session, inputId='feature_picker_gene_module', choices='-spoof-')
    updateSelectizeInput(session=session, inputId='feature_picker_gene_module',
                         choices=feature_picker_options$gene_modules, selected=feature_picker_selected$gene_modules)

    ## hidden picked feature text input
    updateTextInput(session=session, inputId='picked_feature', value='-spoof-')
    updateTextInput(session=session, inputId='picked_feature', value=picked_feature)

    # update the reactive
    previously_picked_feature[[tab]] <- feature_picker_selected})
}
