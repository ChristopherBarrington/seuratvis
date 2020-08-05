#'
#' 
seurat_object_options.ui <- function(id, seurat) {
  pickerInput_defaults <- list(choices=NULL, selected=NULL, multiple=FALSE, inline=FALSE, width=NULL)
  textInput_defaults <- list(value='value', width=NULL, placeholder='placeholder')

  list(inputId=NS(id, 'n_features_picker'), label='Features', options=list(title='Features', size=5)) %>% modifyList(x=pickerInput_defaults) %>% do.call(what=pickerInput) -> n_features_picker
  list(inputId=NS(id, 'n_umi_picker'), label='UMI', options=list(title='UMI', size=5)) %>% modifyList(x=pickerInput_defaults) %>% do.call(what=pickerInput) -> n_umi_picker
  list(inputId=NS(id, 'proportion_mt_picker'), label='Mitochondrial proportion', options=list(title='Mitochondrial proportion', size=5)) %>% modifyList(x=pickerInput_defaults) %>% do.call(what=pickerInput) -> proportion_mt_picker
  list(inputId=NS(id, 'gene_modules_regex_text'), value='^GeneModule-', label='Gene modules regex', placeholder='regex') %>% modifyList(x=textInput_defaults) %>% do.call(what=textInput) -> gene_modules_picker

  # define the biomaRt options
  biomaRt::listEnsemblArchives() %>% filter(str_detect(url, 'archive')) %>% transmute(label=sprintf(fmt='%s [%s]', name, date), value=url) %>% deframe() -> mart_urls
  list(inputId=NS(id, 'mart_url_picker'), label='Ensembl version', options=list(size=5), choices=mart_urls, selected='http://jul2018.archive.ensembl.org') %>%
    modifyList(x=pickerInput_defaults) %>%
    do.call(what=pickerInput) %>%
    conditionalPanel(condition=sprintf('input["%s"]==false', NS(id, 'mart_included'))) -> mart_url_picker

  list(inputId=NS(id, 'mart_species_picker'), label='Species', options=list(title='Species', size=5), choices=c('hsapiens','mmusculus')) %>%
    modifyList(x=pickerInput_defaults) %>%
    do.call(what=pickerInput) %>%
    conditionalPanel(condition=sprintf('input["%s"]==false', NS(id, 'mart_included'))) -> mart_species_picker

  switchInput(inputId=NS(id, 'mart_included'), value=TRUE) %>% hidden() -> mart_included_switch

  # return ui element(s)
  tagList(tags$label('Configure metadata columns'),
          n_features_picker,
          n_umi_picker,
          proportion_mt_picker,
          gene_modules_picker,
          mart_url_picker,
          mart_species_picker,
          mart_included_switch)
}

#'
#' 
seurat_object_choices.ui <- function(id, available_seurats) {
  prettyRadioButtons(inputId=NS(id,'picker'), label='Select a Seurat object',
                     choiceNames=available_seurats$choiceName, choiceValues=available_seurats$choiceValue, selected='',
                     icon=icon('check'), bigger=TRUE, animation='jelly')
}

#'
#' @import purrr
#' @import biomaRt
#' @import gtools
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
    DefaultAssay(s) <- 'RNA'

    # save the Seurat object
    seurat$n_features_updated <- NULL
    seurat$n_umi_updated <- NULL
    seurat$proportion_mt_updated <- NULL

    seurat$object <- s
    seurat$project <- Project(s)
    seurat$formatted_project_name <- seurat$project %>% reformat_project_name()
    seurat$metadata <- s@meta.data
    seurat$n_cells <- nrow(seurat$metadata)
    seurat$features_in_assays <- list()
    seurat$reductions <- Reductions(s)
    seurat$assays <- Assays(s)
    seurat$cluster_resolutions <- c('seurat_clusters', str_subset(colnames(s@meta.data), '_snn_res.'))
    seurat$all_idents <- {resolutions <- c('seurat_clusters', str_subset(colnames(s@meta.data), '_snn_res.')) ; select_at(s@meta.data, vars(all_of(resolutions))) %>% lapply(function(x) {if(is.factor(x)) {levels(x)} else {unique(x)}}) %>% lapply(gtools::mixedsort)}

    # initialise the provenance of the seurat
    seurat$provenance <- s@misc$provenance
    seurat$provenance_missing <- is.null(seurat$provenance)
    if(seurat$provenance_missing) seurat$provenance <- list('aperture'='there is no cake')

    # save biomaRt to reactive, if available
    mart_included <- !is.null(s@misc$mart)
    seurat$mart <- NULL
    if(mart_included)
      seurat$mart <- s@misc$mart
    updateSwitchInput(session=session, inputId='mart_included', value=mart_included)

    # FindMarkers results
    if(!is.null(s@misc$FindMarkersResults$wilcox)) {
      s@misc$FindMarkersResults %>%
        pluck('wilcox') %>%
        mutate(p_adj_group={p_val_adj %>% cut(breaks=c(0, 0.1/100, 1/100, 5/100, 10/100, 100/100), labels=c('<0.1%','<1%','<5%','<10%','NS'), include.lowest=TRUE, right=TRUE)}) %>%
        dplyr::select(cluster_set, ident.1, gene, pct.1, pct.2, avg_logFC, p_adj_group) %>%
        rename(`Cluster set`='cluster_set', `Cluster ID`='ident.1', `Gene`='gene', `Cluster detection`='pct.1', `Map detection`='pct.2', `Avg. logFC`='avg_logFC', `Adj. P`='p_adj_group') -> tidied_results

      ## update the reactive(s)
      seurat$FindMarkersResults$vars <- c('Cluster set', 'Cluster ID', 'Gene', 'Adj. P', 'Avg. logFC', 'Cluster detection', 'Map detection')
      seurat$FindMarkersResults$table <- tidied_results
    } else {
      # seurat$FindMarkersResults$vars <- c('Cluster set', 'Cluster ID', 'Gene', 'Adj. P', 'Avg. logFC', 'Cluster detection', 'Map detection')
      c('Cluster set', 'Cluster ID', 'Gene', 'Adj. P', 'Avg. logFC', 'Cluster detection', 'Map detection') %T>%
        (function(x) seurat$FindMarkersResults$vars <- x) %>%
        purrr::set_names() %>%
        as.list() %>%
        as.data.frame() %>%
        filter(NA) -> seurat$FindMarkersResults$table
    }

    # update ui elements
    ## initialise the metadata selector variables to NULL
    seurat$n_features_variable <- NULL
    seurat$n_umi_variable <- NULL
    seurat$proportion_mt_variable <- NULL

    ## get the numeric metadata variables
    sapply(seurat$metadata, is.numeric) %>% subset(x=colnames(seurat$metadata)) -> numeric_choices
    sapply(seurat$metadata, is.character) %>% subset(x=colnames(seurat$metadata)) -> character_choices

    ## guess a default choice
    n_features_picker_default <- preferred_choice(x=numeric_choices, preferences=c('nFeature_RNA','nFeature_SCT'))
    n_umi_picker_default <- preferred_choice(x=numeric_choices, preferences=c('nCount_RNA','nCount_SCT'))
    proportion_mt_picker_default <- preferred_choice(x=numeric_choices, preferences=c('percent.mt', 'percent_mt', 'prop.mt', 'prop_mt'))

    seurat$n_features_variable <- n_features_picker_default
    seurat$n_umi_variable <- n_umi_picker_default
    seurat$proportion_mt_variable <- proportion_mt_picker_default

    ## define the choices and default in the input ui elements
    updatePickerInput(session=session, inputId='n_features_picker', choices=numeric_choices, selected=n_features_picker_default)
    updatePickerInput(session=session, inputId='n_umi_picker', choices=numeric_choices, selected=n_umi_picker_default)
    updatePickerInput(session=session, inputId='proportion_mt_picker', choices=numeric_choices, selected=proportion_mt_picker_default)})

  # update seurat reactive when a options are selected or object is changed
  ## pick out number of features per cell
  observe(label='process_seurat/n_features_picker', x={
    req(input$n_features_picker)
    req(seurat$metadata)

    if(!is.element(el=input$n_features_picker, set=names(seurat$metadata)))
      return(NULL)

    seurat$n_features_variable <- input$n_features_picker
    seurat$n_features_values <- select(seurat$metadata, input$n_features_picker) %>% unlist(use.names=FALSE)
    seurat$n_features_values_min <- min(seurat$n_features_values)
    seurat$n_features_values_max <- max(seurat$n_features_values)
    seurat$n_features_values_mean <- mean(seurat$n_features_values)
    seurat$n_features_values_median <- median(seurat$n_features_values)

    seurat$n_features_updated <- rnorm(1)})

  ## pick out total umi per cell
  observe(label='process_seurat/n_umi_picker', x={
    req(input$n_umi_picker)
    req(seurat$metadata)

    if(!is.element(el=input$n_umi_picker, set=names(seurat$metadata)))
      return(NULL)

    seurat$n_umi_variable <- input$n_umi_picker
    seurat$n_umi_values <- select(seurat$metadata, input$n_umi_picker) %>% unlist(use.names=FALSE)
    seurat$n_umi_sum <- sum(seurat$n_umi_values)
    seurat$n_umi_values_min <- min(seurat$n_umi_values)
    seurat$n_umi_values_max <- max(seurat$n_umi_values)
    seurat$n_umi_values_mean <- mean(seurat$n_umi_values)
    seurat$n_umi_values_median <- median(seurat$n_umi_values)

    seurat$n_umi_updated <- rnorm(1)})

  ## pick out proportion of mitochondria reads per cell
  observe(label='process_seurat/proportion_mt_picker', x={
    req(input$proportion_mt_picker)
    req(seurat$metadata)

    if(!is.element(el=input$proportion_mt_picker, set=names(seurat$metadata)))
      return(NULL)

    seurat$proportion_mt_variable <- input$proportion_mt_picker
    seurat$proportion_mt_values <- select(seurat$metadata, input$proportion_mt_picker) %>% unlist(use.names=FALSE)
    seurat$proportion_mt_values_min <- min(seurat$proportion_mt_values)
    seurat$proportion_mt_values_max <- max(seurat$proportion_mt_values)
    seurat$proportion_mt_values_mean <- mean(seurat$proportion_mt_values)
    seurat$proportion_mt_values_median <- median(seurat$proportion_mt_values)

    seurat$proportion_mt_updated <- rnorm(1)})

  ## pick out proportion of mitochondria reads per cell
  observe(label='process_seurat/gene_modules_regex', x={
    req(input$gene_modules_regex_text)
    req(seurat$metadata)

    gm_regex <- input$gene_modules_regex_text

    seurat$gene_modules <- ''

    if(!is.null(seurat$object@misc$gene_modules))
      seurat$gene_modules <- seurat$object@misc$gene_modules %>% purrr::set_names(str_remove, pattern=regex(gm_regex))
    seurat$gene_module_scores <- dplyr::select(seurat$metadata, matches(gm_regex)) %>% purrr::set_names(str_remove, pattern=regex(gm_regex))
    seurat$gene_modules_regex <- gm_regex})

  ## configure the biomaRt object
  observe(label='process_seurat/update_mart_species', x={
    req(input$mart_url_picker)

    if(!input$mart_included)
      biomaRt::useMart(biomart='ensembl', host=input$mart_url_picker) %>%
        biomaRt::listDatasets() %>%
        pluck('dataset') %>%
        str_remove('_gene_ensembl') %>%
        updatePickerInput(session=session, inputId='mart_species_picker', label=NULL, selected=NULL)})

  observe(label='process_seurat/configure_biomart', x={
    req(input$mart_url_picker)
    req(input$mart_species_picker)

    if(!input$mart_included)
      biomaRt::useMart(biomart='ensembl', host=input$mart_url_picker) %>%
        biomaRt::useDataset(dataset=sprintf('%s_gene_ensembl', input$mart_species_picker)) -> seurat$mart})

  # return the reactive
  seurat
}
