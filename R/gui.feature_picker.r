#' 
#' 
feature_picker.ui <- function(id, seurat, label='Feature selection', selected='features', include_feature_type=TRUE, include_values_range=TRUE,
                              choices=list(`Features`='features', `Metadata`='metadata', `Gene modules`='gene_modules'),
                              features_opts=list(), metadata_opts=list(), gene_modules_opts=list(),
                              features_regex='.*', metadata_regex='.*', gene_modules_regex='.*',
                              metadata_filter=function(x) x) {
  ns <- NS(id)

  # get the possible features and values
  ## get names of features and metadata
  list(features=rownames(seurat$object),
       metadata=seurat$metadata %>% metadata_filter() %>% colnames(),
       gene_modules=colnames(seurat$gene_module_scores)) -> feature_picker_options

  ## filter the list for non-empty sets
  feature_picker_options <- feature_picker_options[sapply(feature_picker_options, length)>0]
  
  ## only use choices with non-empty option sets
  choices <- choices[unlist(choices) %in% names(feature_picker_options)]

  ## filter the options using the regex
  feature_picker_options$features %<>% str_subset(pattern=regex(pattern=features_regex, ignore_case=TRUE))
  feature_picker_options$metadata %<>% str_subset(pattern=regex(pattern=metadata_regex, ignore_case=TRUE))

  ## pick a random feature and metadata and gene module column
  feature_picker_options %>%
    lapply(sample, size=1) -> feature_picker_selected

  # pick a feature to display: features or metadata or gene_module
  picked_feature <- feature_picker_selected$features

  # make ui elements
  ## feature names autocomplete box
  autocomplete_input(id=ns(id='feature_picker_feature_names'), label=NULL, placeholder='Feature',
                     options=feature_picker_options$features, value=feature_picker_selected$features) %>%
    conditionalPanel(condition=sprintf('input["%s"]=="features"', ns(id='feature_type'))) -> feature_names_picker_conditional

  ## metadata names drop down box
  # selectizeInput(inputId=ns(id='feature_picker_metadata'), label=NULL,
  #                choices=feature_picker_options$metadata, selected=feature_picker_selected$metadata, multiple=FALSE) %>%
  #   conditionalPanel(condition=sprintf('input["%s"]=="metadata"', ns(id='feature_type'))) -> metadata_picker_conditional

  list(inputId=ns(id='feature_picker_metadata'), label=NULL,
       choices=feature_picker_options$metadata, selected=feature_picker_selected$metadata, multiple=FALSE) %>%
    modifyList(val=metadata_opts) %>%
    do.call(what=selectizeInput) %>%
    conditionalPanel(condition=sprintf('input["%s"]=="metadata"', ns(id='feature_type'))) -> metadata_picker_conditional

  ## gene modules drop down box
  list(inputId=ns(id='feature_picker_gene_module'), label=NULL,
       choices=feature_picker_options$gene_modules, selected=feature_picker_options$gene_modules, multiple=FALSE,
       options=list(`actions-box`=TRUE, header='Gene module(s) selection', title='Gene module selection',
                    `selected-text-format`='count', `count-selected-text`='{0} module(s)')) %>%
    modifyList(val=gene_modules_opts) %>%
    do.call(what=pickerInput) %>%
    conditionalPanel(condition=sprintf('input["%s"]=="gene_modules"', ns(id='feature_type'))) -> gene_module_picker_conditional

  ## slider to limit colour range
  sliderInput(inputId=ns(id='value_range'), label='Colour range limits',
              min=0, max=1, step=0.1, value=c(-Inf,Inf)) -> value_range

  ## checkbox for feature type
  prettyRadioButtons(inputId=ns(id='feature_type'), status='primary', label=label, 
                     choices=choices, selected=selected,
                     icon=icon('check'), bigger=TRUE, animation='jelly') -> feature_type_picker
  if(!include_feature_type)
    feature_type_picker %<>% hidden()

  ## hidden text box to serve app
  textInput(inputId=ns('picked_feature'), label='picked feature', value=picked_feature) %>% hidden() -> picked_feature_text_input

  # return ui element(s)
  tagList(feature_type_picker,
          feature_names_picker_conditional,
          metadata_picker_conditional,
          gene_module_picker_conditional,
          if(include_values_range) value_range,
          picked_feature_text_input)
}

#' 
#' 
feature_picker.server <- function(input, output, session, seurat, features_regex='.*', metadata_regex='.*', ...) {
  # previously_picked_feature <- reactiveValues()
  picked_feature <- reactiveValues()

  # react to the feature selection
  ## if a feature is selected, copy it to the reactive
  observeEvent(eventExpr=input$feature_picker_feature_names, handlerExpr={
    # make sure these elements are defined
    req(input$feature_picker_feature_names)
   
    # update hidden ui element
    if(input$feature_type=='features')
      picked_feature$name <- input$feature_picker_feature_names})
      # updateTextInput(session=session, inputId='picked_feature', value=input$feature_picker_feature_names)})

  ## if a metadata column is selected, copy it to the reactive
  observeEvent(eventExpr=input$feature_picker_metadata, handlerExpr={
    # make sure these elements are defined
    req(input$feature_picker_metadata)

    # update hidden ui element
    if(input$feature_type=='metadata')
      picked_feature$name <- input$feature_picker_metadata})
      # updateTextInput(session=session, inputId='picked_feature', value=input$feature_picker_metadata)})

  ## if a gene module column is selected, copy it to the reactive
  observeEvent(eventExpr=input$feature_picker_gene_module, handlerExpr={
    # make sure these elements are defined
    req(input$feature_picker_gene_module)

    # update hidden ui element
    if(input$feature_type=='gene_modules')
      picked_feature$name <- input$feature_picker_gene_module})
      # updateTextInput(session=session, inputId='picked_feature', value=input$feature_picker_gene_module)})

  ## update the hidden ui element when a feature type is selected
  observeEvent(eventExpr=input$feature_type, handlerExpr={
    # pick the feature to revert to
    input_name <- switch(input$feature_type, features='feature_picker_feature_names', metadata='feature_picker_metadata', gene_modules='feature_picker_gene_module')

    # update hidden ui element
    picked_feature$name <- input[[input_name]]})
    # updateTextInput(session=session, inputId='picked_feature', value=input[[input_name]])})

  ## use the selected feature (it may be a feature or metadata)
  # observeEvent(eventExpr=input$picked_feature, handlerExpr={
  observe(label='feature_picker/fetch', x={
    # make sure these elements are defined
    req(seurat$object)
    req(input$feature_type)
    req(picked_feature$name)

    # create variables for shorthand
    picked <- picked_feature$name

    # get the values for the selected feature(s) from the loaded Seurat
    if(input$feature_type=='gene_modules') {
      picked %<>% str_split(pattern=',') %>% unlist()
      picked_feature_values <- select_at(seurat$gene_module_scores, vars(picked))
    } else {
      picked_feature_values <- FetchData(object=seurat$object, vars=picked)
    }

    if(length(picked)==1) {
      picked_feature_values %<>% set_names('value')

      # update the ui element(s)
      ## slider to limit colour range
      min_value <- 0
      max_value <- 1
      if(!is.null(picked_feature_values$value) && class(picked_feature_values$value)=='numeric') {
        min_value <- min(picked_feature_values$value) %>% subtract(0.05) %>% round(digits=1)
        max_value <- max(picked_feature_values$value) %>% add(0.05) %>% round(digits=1)
      }

      updateSliderInput(session=session, inputId='value_range',
                        min=min_value, max=max_value, value=c(-Inf,Inf))
    } else {
      updateSliderInput(session=session, inputId='value_range',
                        min=0, max=0, value=c(-Inf,Inf))
    }

    # save feature information in the reactive
    picked_feature$values <- picked_feature_values

    # invalidate the reactive when both data and slider are updated
    picked_feature$refreshed <- rnorm(1)})

  # invalidate the reactive value when slider is changed but not after initialisation
  observeEvent(eventExpr=input$value_range, ignoreInit=TRUE, handlerExpr={
    picked_feature$refreshed <- rnorm(1)
    picked_feature$values_range <- input$value_range})

  # reset the reactive when the seurat is (re)loaded
  observeEvent(eventExpr=seurat$object, handlerExpr={
    for(i in names(picked_feature))
      picked_feature[[i]] <- NULL})

  # return the reactiveValues list
  return(picked_feature)
}
