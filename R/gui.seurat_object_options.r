#'
#' 
seurat_object_options.ui <- function(id, seurat) {
  pickerInput_defaults <- list(choices=NULL, selected=NULL, multiple=FALSE, inline=FALSE, width=NULL)
  list(inputId=NS(id, 'n_features_picker'), label='Features', options=list(title='Features')) %>% modifyList(x=pickerInput_defaults) %>% do.call(what=pickerInput) -> n_features_picker
  list(inputId=NS(id, 'n_umi_picker'), label='UMI', options=list(title='UMI')) %>% modifyList(x=pickerInput_defaults) %>% do.call(what=pickerInput) -> n_umi_picker
  list(inputId=NS(id, 'proportion_mt_picker'), label='Mitochondrial proportion', options=list(title='Mitochondrial proportion')) %>% modifyList(x=pickerInput_defaults) %>% do.call(what=pickerInput) -> proportion_mt_picker

  # return ui element(s)
  tagList(tags$label('Configure metadata columns'),
          n_features_picker,
          n_umi_picker,
          proportion_mt_picker)
}

#'
#' 
seurat_object_options.server <- function(input, output, session, server_input, server_session, seurat) {
  
  # update the ui elements when an object is loaded
  observeEvent(eventExpr=seurat$metadata, label='seurat_object_options/object', handlerExpr={
    # update ui elements
    ## get the numeric metadata variables
    sapply(seurat$metadata, is.numeric) %>% subset(x=colnames(seurat$metadata)) -> choices

    ## guess a default choice
    n_features_picker_default <- preferred_choice(x=choices, preferences=c('nFeature_RNA','nFeature_SCT'))
    n_umi_picker_default <- preferred_choice(x=choices, preferences=c('nCount_RNA','nCount_SCT'))
    proportion_mt_picker_default <- preferred_choice(x=choices, preferences=c('percent.mt', 'percent_mt', 'prop.mt', 'prop_mt'))

    ## define the choices and default in the input ui elements
    updatePickerInput(session=server_session, inputId=NS('process_seurat', 'n_features_picker'), choices=choices, selected=n_features_picker_default)
    updatePickerInput(session=server_session, inputId=NS('process_seurat', 'n_umi_picker'), choices=choices, selected=n_umi_picker_default)
    updatePickerInput(session=server_session, inputId=NS('process_seurat', 'proportion_mt_picker'), choices=choices, selected=proportion_mt_picker_default)})

  # update seurat reactive when a column is selected
  observeEvent(eventExpr=server_input[[NS('process_seurat', 'n_features_picker')]], label='', handlerExpr={seurat$n_features_variable <- server_input[[NS('process_seurat', 'n_features_picker')]]})
  observeEvent(eventExpr=server_input[[NS('process_seurat', 'n_umi_picker')]], label='', handlerExpr={seurat$n_umi_variable <- server_input[[NS('process_seurat', 'n_umi_picker')]]})
  observeEvent(eventExpr=server_input[[NS('process_seurat', 'proportion_mt_picker')]], label='', handlerExpr={seurat$proportion_mt_variable <- server_input[[NS('process_seurat', 'proportion_mt_picker')]]})
}
