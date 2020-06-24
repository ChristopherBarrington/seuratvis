#'
#' 
cluster_classification.tab <- function() {
  bquote({
    menuItem(text='Cluster classification', icon=icon('object-group'), startExpanded=TRUE,
             menuSubItem(text='FindMarkers results', tabName='findmarkers_results_tab'),
             menuItem(text='Gene modules', startExpanded=TRUE,
                      menuSubItem(text='Score in clusters', tabName='gene_module_score_in_clusters_tab'),
                      menuSubItem(text='Scores in a cluster', tabName='gene_module_scores_in_a_cluster_tab')) %>%
              modify_stop_propagation()) -> menu_item

    list(tabItem(tabName='findmarkers_results_tab',
                 h1('TITLE'),
                 fluidRow(project_name_text_box.ui(id=NS('findmarkers_results_tab', 'project_name'), width=12)),
                 fluidRow(boxPlus(title='Results table', closable=FALSE, width=12, status='primary'))),
         tabItem(tabName='gene_module_score_in_clusters_tab',
                 h1('Gene module score in clusters'),
                 fluidRow(project_name_text_box.ui(id=NS('gene_module_score_in_clusters_tab', 'project_name'), width=12)),
                 column(width=7,
                        boxPlus(title='Gene module score', closable=FALSE, width=12, height='75vh', status='primary', 
                                feature_ridges.plot(id=NS('gene_module_score_in_clusters_tab', 'scores_plot')))),
                 column(width=5,
                        boxPlus(title='Map', closable=FALSE, width=12, height='35vh', status='primary',
                                dimension_reduction.plot(id=NS('gene_module_score_in_clusters_tab', 'map'))),
                        boxPlus(title='Genes', closable=FALSE, width=12, height='35vh', status='primary',
                                genes_in_modules.table(id=NS('gene_module_score_in_clusters_tab', 'modules_table'))))),
         tabItem(tabName='gene_module_scores_in_a_cluster_tab',
                 h1('Gene modules'),
                 fluidRow(project_name_text_box.ui(id=NS('gene_module_scores_in_a_cluster_tab', 'project_name'), width=12)),
                 fluidRow(column(width=12, feature_ridges.plot(id=NS('gene_module_scores_in_a_cluster_tab', 'scores_plot')))))) -> content

    menus %<>% append(list(menu_item))
    contents %<>% append(content)})
}

#'
#' 
findmarkers_results_tab.server <- function(input, output, session, server_input, server_output, server_session, seurat) {
  # build the sidebar ui
  observeEvent(eventExpr=server_input$left_sidebar, handlerExpr={
    tab <- 'findmarkers_results_tab'
    if(server_input$left_sidebar==tab) {
      tab %<>% str_c('-')
      renderUI({tagList()})  -> server_output$right_sidebar.data_opts
      renderUI({tagList()}) -> server_output$right_sidebar.plotting_opts}})

  # call the modules for this tab
  callModule(module=project_name_text_box.server, id='project_name', seurat=seurat)
}

#'
#' 
gene_module_score_in_clusters_tab.server <- function(input, output, session, server_input, server_output, server_session, seurat) {
  # build the sidebar ui
  observeEvent(eventExpr=server_input$left_sidebar, handlerExpr={
    print(server_input$left_sidebar)
    tab <- 'gene_module_score_in_clusters_tab'
    if(server_input$left_sidebar==tab) {
      tab %<>% str_c('-')
      renderUI({tagList(cluster_picker.ui(id=tab, seurat=seurat, resolution=TRUE, picker=TRUE, label_switch=FALSE),
                        feature_picker.ui(id=tab, seurat=seurat, choices=list(`Gene modules`='gene_modules'), selected='gene_modules', gene_modules_opts=list(multiple=FALSE), include_feature_type=FALSE, include_values_range=FALSE),
                        dimension_reduction.ui(id=tab, seurat=seurat))})  -> server_output$right_sidebar.data_opts
      renderUI({tagList()}) -> server_output$right_sidebar.plotting_opts}})

  # call the modules for this tab
  dimension_reduction <- callModule(module=dimension_reduction.server, id='', seurat=seurat, regex='.*')
  cluster_resolution <- callModule(module=cluster_picker.server, id='', seurat=seurat)
  feature_picker <- callModule(module=feature_picker.server, id='', seurat=seurat)
  colour_picker <- list(low='linen', mid='white', high='darkviolet', background=rgb(255, 255, 255, 255, max=255))

  callModule(module=project_name_text_box.server, id='project_name', seurat=seurat)
  callModule(module=feature_ridge_by_idents.server, id='scores_plot', picked_feature=feature_picker, picked_clusters=cluster_resolution)
  callModule(module=dimension_reduction.show_selected_clusters.server, id='map', dimension_reduction=dimension_reduction, point_size=list(size=0.6), cluster_resolution=cluster_resolution, picked_colours=colour_picker)
  callModule(module=genes_in_modules.server, id='modules_table', seurat=seurat, picked_feature=feature_picker)
}

