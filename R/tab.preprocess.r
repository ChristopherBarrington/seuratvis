
preprocessing.tab <- function() {
  bquote({
    menuItem(text='Preprocessing', icon=icon('toolbox'), startExpanded=TRUE,
             menuSubItem(text='Cell filtering', tabName='cell_filtering_tab'),
             menuSubItem(text='Dimensionality', tabName='dimensionality_tab')) -> menu_item

    list(tabItem(tabName='cell_filtering_tab',
                 h1('Identify cells that can be removed with quality filters'),
                 fluidRow(project_name_text_box.ui(id=NS('cell_filtering_tab', 'project_name'), width=12)),
                 fluidRow(boxPlus(title='UMI per cell', closable=FALSE, width=12, status='primary'
                                  # boxplot.plot(id='cell_filtering_tab-umi_per_cell'),
                                  # knee.plot(id='cell_filtering_tab-umi_per_cell'),
                                  # density_line.plot(id='cell_filtering_tab-umi_per_cell')
                                  )),
                 fluidRow(boxPlus(title='Features per cell', closable=FALSE, width=12, status='primary'
                                  # boxplot.plot(id='cell_filtering_tab-picked_feature'),
                                  # knee.plot(id='cell_filtering_tab-picked_feature'),
                                  # density_line.plot(id='cell_filtering_tab-picked_feature')
                                  )),
                 fluidRow(boxPlus(title='Mitochondrial UMI', closable=FALSE, width=12, status='primary'
                                  # boxplot.plot(id='cell_filtering_tab-mt_umi'),
                                  # knee.plot(id='cell_filtering_tab-mt_umi'),
                                  # density_line.plot(id='cell_filtering_tab-mt_umi')
                                  ))),
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
      addClass(selector='body', class='control-sidebar-open')
      showTab(inputId='right_sidebar_tabset', target='data_opts', select=TRUE, session=server_session)
      renderUI({tagList(filter_umi_per_cell.ui(id=tab, seurat=seurat),
                        filter_features_per_cell.ui(id=tab, seurat=seurat, resolution=TRUE, label_switch=TRUE),
                        filter_mt_umi.ui(id=tab, seurat=seurat),
                        show_filtering_paramters.ui())})  -> server_output$right_sidebar.data_opts
      renderUI({tagList()}) -> server_output$right_sidebar.plotting_opts}})

  # call the modules for this tab
  callModule(module=project_name_text_box.server, id='project_name', seurat=seurat)
}

dimensionality_tab.server <- function(input, output, session, server_input, server_output, server_session, seurat) {
  # build the sidebar ui
  observeEvent(eventExpr=server_input$left_sidebar, handlerExpr={
    tab <- 'dimensionality_tab'
    if(server_input$left_sidebar==tab) {
      tab %<>% str_c('-')
      addClass(selector='body', class='control-sidebar-open')
      showTab(inputId='right_sidebar_tabset', target='data_opts', select=TRUE, session=server_session)
      renderUI({tagList()})  -> server_output$right_sidebar.data_opts
      renderUI({tagList()}) -> server_output$right_sidebar.plotting_opts}})

  # call the modules for this tab
  callModule(module=project_name_text_box.server, id='project_name', seurat=seurat)
}