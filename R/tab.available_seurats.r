#'
#' 
available_seurats.tab <- function() {
  bquote({

    tab <- 'configuration_tab'
    menuItem(text='Configure', tabName=tab, icon=icon('map'), selected=TRUE) -> menu_item

    tabItem(tabName='configuration_tab',
            h1('Available Seurat objects'),
            h5('These objects have been found in your workspace and can be loaded'),
            fluidRow(boxPlus(title='', closable=FALSE, width=12, status='primary',
                             seurats_in_workspace.table(id=NS(tab, 'seurats_table'))))) -> content

    menus %<>% append(list(menu_item))
    contents %<>% append(list(content))})
}

#'
#' 
available_seurats_tab.server <- function(input, output, session, server_input, server_output, server_session) {
  # build the sidebar ui
  observeEvent(eventExpr=server_input$left_sidebar, handlerExpr={
    tab <- 'configuration_tab'
    if(server_input$left_sidebar==tab) {
      tab %<>% str_c('-')
      renderUI({p('No options')})  -> server_output$right_sidebar.data_opts
      renderUI({p('No options')}) -> server_output$right_sidebar.plotting_opts}})

  # call the modules for this tab
  callModule(seurats_in_workspace.server, id='seurats_table')
}

#'
#' 
process_seurat.server <- function(input, output, session, server_input, server_output, server_session, available_seurats) {
  seurat <- reactiveValues()

  # render the config right sidebar tab when the server starts
  renderUI({tagList(prettyRadioButtons(inputId=NS('process_seurat','picker'), label='Available Seurat objects!',
                                       choiceNames=available_seurats$choiceName, choiceValues=available_seurats$choiceValue, selected='',
                                       icon=icon('check'), bigger=TRUE, animation='jelly'),
                    seurat_object_options.ui(id='process_seurat', seurat=seurat))}) -> server_output$right_sidebar.config_opts

  callModule(module=seurat_object_options.server, id='process_seurat', server_input=server_input, server_session=server_session, seurat=seurat)

  # react when a seurat is selected
  observeEvent(server_input[['process_seurat-picker']], {
    s <- eval(parse(text=server_input[['process_seurat-picker']]))
    seurat$object <- s
    seurat$formatted_project_name <- Project(s) %>% reformat_project_name()
    seurat$metadata <- s@meta.data
    seurat$features_in_assays <- list()
    seurat$reductions <- Reductions(s)
    seurat$assays <- Assays(s)
    seurat$gene_module_scores <- select_at(s@meta.data, vars(starts_with('GeneModule-')))
    seurat$gene_modules <- s@misc$gene_modules
    seurat$cluster_resolutions <- c('seurat_clusters', str_subset(colnames(s@meta.data), '_snn_res.'))
    seurat$all_idents <- {resolutions <- c('seurat_clusters', str_subset(colnames(s@meta.data), '_snn_res.')) ; select_at(s@meta.data, vars(all_of(resolutions))) %>% plyr::llply(levels)}
  })

  # return the reactive
  seurat
}
