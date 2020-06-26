
preprocessing.tab <- function() {
  bquote({
    menuItem(text='Preprocessing', icon=icon('toolbox'), startExpanded=TRUE,
             menuSubItem(text='Cell filtering', tabName='cell_filtering_tab'),
             menuSubItem(text='Dimensionality', tabName='dimensionality_tab')) -> menu_item

    list(tabItem(tabName='cell_filtering_tab',
                 h1('Identify cells that can be removed with quality filters'),
                 fluidRow(project_name_text_box.ui(id=NS('cell_filtering_tab', 'project_name'), width=12)),
                 fluidRow(boxPlus(title='UMI per cell', closable=FALSE, width=12, status='primary',
                                  column(width=2, offset=1, boxplot.plot(id=NS('cell_filtering_tab', 'n_umi'))),
                                  column(width=3, offset=0, knee.plot(id=NS('cell_filtering_tab', 'n_umi'))),
                                  column(width=4, offset=1, density.plot(id=NS('cell_filtering_tab', 'n_umi'))))),
                 fluidRow(boxPlus(title='Features per cell', closable=FALSE, width=12, status='primary',
                                  column(width=2, offset=1, boxplot.plot(id=NS('cell_filtering_tab', 'n_features'))),
                                  column(width=3, offset=0, knee.plot(id=NS('cell_filtering_tab', 'n_features'))),
                                  column(width=4, offset=1, density.plot(id=NS('cell_filtering_tab', 'n_features'))))),
                 fluidRow(boxPlus(title='Mitochondrial UMI', closable=FALSE, width=12, status='primary',
                                  column(width=2, offset=1, boxplot.plot(id=NS('cell_filtering_tab', 'proportion_mt'))),
                                  column(width=3, offset=0, knee.plot(id=NS('cell_filtering_tab', 'proportion_mt'))),
                                  column(width=4, offset=1, density.plot(id=NS('cell_filtering_tab', 'proportion_mt')))))),
         tabItem(tabName='dimensionality_tab',
                 h1('Determine the dimensionality of a dataset'),
                 fluidRow(project_name_text_box.ui(id=NS('dimensionality_tab', 'project_name'), width=12)),
                 fluidRow(boxPlus(title='Elbow plot', closable=FALSE, width=6, status='primary'),
                          boxPlus(title='Top contributing features', closable=FALSE, width=6, status='primary'),
                          boxPlus(title='JackStraw plot', closable=FALSE, width=6, status='primary'),
                          boxPlus(title='JackStraw p-values', closable=FALSE, width=6, status='primary')))) -> content

    menus %<>% append(list(menu_item))
    contents %<>% append(content)})
}

cell_filtering_tab.server <- function(input, output, session, server_input, server_output, server_session, seurat) {
  # build the sidebar ui
  observeEvent(eventExpr=server_input$left_sidebar, handlerExpr={
    tab <- 'cell_filtering_tab'
    if(server_input$left_sidebar==tab) {    
      tab %<>% str_c('-')
      renderUI({tagList(filter_n_umi.ui(id=tab, seurat=seurat),
                        filter_n_features.ui(id=tab, seurat=seurat),
                        filter_proportion_mt.ui(id=tab, seurat=seurat),
                        show_filtering_parameters.ui(id=tab, label='Cell filtering parameters'))})  -> server_output$right_sidebar.data_opts
      renderUI({tagList()}) -> server_output$right_sidebar.plotting_opts}})

  # call the modules for this tab
  filtering_parameters <- callModule(module=dataset_filtering.server, id='filtering', seurat=seurat)

  callModule(module=project_name_text_box.server, id='project_name', seurat=seurat)

  callModule(module=filter_n_umi.server, id='', cell_filtering=filtering_parameters)
  callModule(module=filter_n_features.server, id='', cell_filtering=filtering_parameters)
  callModule(module=filter_proportion_mt.server, id='', cell_filtering=filtering_parameters)

  callModule(module=boxplot_plot.n_features.server, id='n_features', seurat=seurat, cell_filtering=filtering_parameters)
  callModule(module=boxplot_plot.n_umi.server, id='n_umi', seurat=seurat, cell_filtering=filtering_parameters)
  callModule(module=boxplot_plot.proportion_mt.server, id='proportion_mt', seurat=seurat, cell_filtering=filtering_parameters)

  callModule(module=knee_plot.n_features.server, id='n_features', seurat=seurat, cell_filtering=filtering_parameters)
  callModule(module=knee_plot.n_umi.server, id='n_umi', seurat=seurat, cell_filtering=filtering_parameters)
  callModule(module=knee_plot.proportion_mt.server, id='proportion_mt', seurat=seurat, cell_filtering=filtering_parameters)

  callModule(module=density_plot.n_features.server, id='n_features', seurat=seurat, cell_filtering=filtering_parameters)
  callModule(module=density_plot.n_umi.server, id='n_umi', seurat=seurat, cell_filtering=filtering_parameters)
  callModule(module=density_plot.proportion_mt.server, id='proportion_mt', seurat=seurat, cell_filtering=filtering_parameters)

  callModule(module=show_filtering_parameters.server, id='', server_output=server_output, seurat=seurat, cell_filtering=filtering_parameters)
}

dimensionality_tab.server <- function(input, output, session, server_input, server_output, server_session, seurat) {
  # build the sidebar ui
  observeEvent(eventExpr=server_input$left_sidebar, handlerExpr={
    tab <- 'dimensionality_tab'
    if(server_input$left_sidebar==tab) {
      tab %<>% str_c('-')
      renderUI({tagList()})  -> server_output$right_sidebar.data_opts
      renderUI({tagList()}) -> server_output$right_sidebar.plotting_opts}})

  # call the modules for this tab
  callModule(module=project_name_text_box.server, id='project_name', seurat=seurat)
}