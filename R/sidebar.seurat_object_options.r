#'
#' 
seurat_object_options.ui <- function(id, seurat) {
  pickerInput_defaults <- list(choices=NULL, selected=NULL, multiple=FALSE, inline=FALSE, width=NULL)
  textInput_defaults <- list(value='value', width=NULL, placeholder='placeholder')

  list(inputId=NS(id, 'n_features_picker'), label='Features', options=list(title='Features')) %>% modifyList(x=pickerInput_defaults) %>% do.call(what=pickerInput) -> n_features_picker
  list(inputId=NS(id, 'n_umi_picker'), label='UMI', options=list(title='UMI')) %>% modifyList(x=pickerInput_defaults) %>% do.call(what=pickerInput) -> n_umi_picker
  list(inputId=NS(id, 'proportion_mt_picker'), label='Mitochondrial proportion', options=list(title='Mitochondrial proportion')) %>% modifyList(x=pickerInput_defaults) %>% do.call(what=pickerInput) -> proportion_mt_picker
  list(inputId=NS(id, 'gene_modules_regex_text'), value='^GeneModule-', label='Gene modules regex', placeholder='regex') %>% modifyList(x=textInput_defaults) %>% do.call(what=textInput) -> gene_modules_picker

  # return ui element(s)
  tagList(tags$label('Configure metadata columns'),
          n_features_picker,
          n_umi_picker,
          proportion_mt_picker,
          gene_modules_picker)
}

#'
#' 
seurat_object_choices.ui <- function(id, available_seurats) {
  prettyRadioButtons(inputId=NS(id,'picker'), label='Select a Seurat object',
                     choiceNames=available_seurats$choiceName, choiceValues=available_seurats$choiceValue, selected='',
                     icon=icon('check'), bigger=TRUE, animation='jelly')
}

#'
#' 
process_seurat.server <- function(input, output, session, server_input, server_output, server_session, available_seurats) {
  seurat <- reactiveValues()

  # render the config right sidebar tab when the server starts
  renderUI({tagList(seurat_object_choices.ui(id='process_seurat', available_seurats=available_seurats),
                    seurat_object_options.ui(id='process_seurat', seurat=seurat))}) -> server_output$right_sidebar.config_opts

  # callModule(module=seurat_object_options.server, id='process_seurat', server_input=server_input, server_session=server_session, seurat=seurat)

  # react when a seurat is selected
  observeEvent(eventExpr=input$picker, label='process_seurat/picker', handlerExpr={
    s <- eval(parse(text=input$picker))

    # check that there are clusters, if not add a fake one
    if(is.null(s@meta.data$seurat_clusters))
      s@meta.data$seurat_clusters <- 0

    # ensure we are using RNA assay
    DefaultAssay(s) <- selected_assay

    seurat$object <- s
    seurat$formatted_project_name <- Project(s) %>% reformat_project_name()
    seurat$metadata <- s@meta.data
    seurat$features_in_assays <- list()
    seurat$reductions <- Reductions(s)
    seurat$assays <- Assays(s)
    seurat$gene_modules <- s@misc$gene_modules
    seurat$cluster_resolutions <- c('seurat_clusters', str_subset(colnames(s@meta.data), '_snn_res.'))
    seurat$all_idents <- {resolutions <- c('seurat_clusters', str_subset(colnames(s@meta.data), '_snn_res.')) ; select_at(s@meta.data, vars(all_of(resolutions))) %>% plyr::llply(levels)}

    # update ui elements
    ## get the numeric metadata variables
    sapply(seurat$metadata, is.numeric) %>% subset(x=colnames(seurat$metadata)) -> numeric_choices
    sapply(seurat$metadata, is.character) %>% subset(x=colnames(seurat$metadata)) -> character_choices

    ## guess a default choice
    n_features_picker_default <- preferred_choice(x=numeric_choices, preferences=c('nFeature_RNA','nFeature_SCT'))
    n_umi_picker_default <- preferred_choice(x=numeric_choices, preferences=c('nCount_RNA','nCount_SCT'))
    proportion_mt_picker_default <- preferred_choice(x=numeric_choices, preferences=c('percent.mt', 'percent_mt', 'prop.mt', 'prop_mt'))

    ## define the choices and default in the input ui elements
    updatePickerInput(session=session, inputId='n_features_picker', choices=numeric_choices, selected=n_features_picker_default)
    updatePickerInput(session=session, inputId='n_umi_picker', choices=numeric_choices, selected=n_umi_picker_default)
    updatePickerInput(session=session, inputId='proportion_mt_picker', choices=numeric_choices, selected=proportion_mt_picker_default)})

  # update seurat reactive when a options are selected or object is changed
  ## pick out number of features per cell
  observe(label='process_seurat/n_features_picker', x={
    req(input$n_features_picker)
    req(seurat$metadata)

    seurat$n_features_variable <- input$n_features_picker
    seurat$n_features_values <- select(seurat$metadata, input$n_features_picker) %>% unlist(use.names=FALSE)
    seurat$n_features_min <- min(seurat$n_features_values)
    seurat$n_features_values_max <- max(seurat$n_features_values)
    seurat$n_features_mean <- mean(seurat$n_features_values)
    seurat$n_features_median <- median(seurat$n_features_values)})

  ## pick out total umi per cell
  observe(label='process_seurat/n_umi_picker', x={
    req(input$n_umi_picker)
    req(seurat$metadata)

    seurat$n_umi_variable <- input$n_umi_picker
    seurat$n_umi_values <- select(seurat$metadata, input$n_umi_picker) %>% unlist(use.names=FALSE)
    seurat$n_umi_min <- min(seurat$n_umi_values)
    seurat$n_umi_values_max <- max(seurat$n_umi_values)
    seurat$n_umi_mean <- mean(seurat$n_umi_values)
    seurat$n_umi_median <- median(seurat$n_umi_values)})

  ## pick out proportion of mitochondria reads per cell
  observe(label='process_seurat/proportion_mt_picker', x={
    req(input$proportion_mt_picker)
    req(seurat$metadata)

    seurat$proportion_mt_variable <- input$proportion_mt_picker
    seurat$proportion_mt_values <- select(seurat$metadata, input$proportion_mt_picker) %>% unlist(use.names=FALSE)
    seurat$proportion_mt_min <- min(seurat$proportion_mt_values)
    seurat$proportion_mt_values_max <- max(seurat$proportion_mt_values)
    seurat$proportion_mt_mean <- mean(seurat$proportion_mt_values)
    seurat$proportion_mt_median <- median(seurat$proportion_mt_values)})

  ## pick out proportion of mitochondria reads per cell
  observe(label='process_seurat/proportion_mt_picker', x={
    req(input$gene_modules_regex_text)
    req(seurat$metadata)

    seurat$gene_modules_regex <- input$gene_modules_regex_text
    seurat$gene_module_scores <- select(seurat$metadata, matches(input$gene_modules_regex_text))})

  # return the reactive
  seurat
}
