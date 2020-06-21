
highlight_features.tab <- function() {
  bquote({
    menuItem(text='Highlight features', icon=icon('chart-line'),
             menuSubItem(text='Highlight feature', tabName='highlight_feature_tab'),
             menuSubItem(text='Highlight feature and cluster', tabName='highlight_feature_and_clusters_tab')) -> menu_item

    list(tabItem(tabName='highlight_feature_tab',
                 h1('Highlight a feature on a map'),
                 fluidRow(boxPlus(title='All clusters in map', closable=FALSE, width=4, dimension_reduction.plot(id='highlight_feature_tab-all_clusters')),
                          boxPlus(title='Selected feature', closable=FALSE, width=4, dimension_reduction.plot(id='highlight_feature_tab-picked_feature')),
                          boxPlus(title='Feature values in clusters', closable=FALSE, width=4, feature_value_per_cluster.plot(id='highlight_feature_tab-feature_value_per_cluster')))),
         tabItem(tabName='highlight_feature_and_clusters_tab',
                 h1('Highlight a feature and clusters'),
                 fluidRow(boxPlus(title='All clusters in map', closable=FALSE, width=4, dimension_reduction.plot(id='highlight_feature_and_clusters_tab-all_clusters')),
                          boxPlus(title='Selected feature', closable=FALSE, width=4, dimension_reduction.plot(id='highlight_feature_and_clusters_tab-picked_feature')),
                          boxPlus(title='Feature values in clusters', closable=FALSE, width=4, feature_value_per_cluster.plot(id='highlight_feature_and_clusters_tab-feature_value_per_cluster'))))) -> content

      menus %<>% append(list(menu_item))
      contents %<>% append(content)})
}

highlight_feature_tab.server <- function(input, output, session, server_input, server_output, server_session, seurat) {
  # build the sidebar ui
  tab <- 'highlight_feature_tab-'
  observeEvent(eventExpr=server_input$left_sidebar, handlerExpr={
    if(server_input$left_sidebar=='highlight_feature_tab') {
      shinyjs::addClass(selector='body', class='control-sidebar-open')
      renderUI({tagList(dimension_reduction.ui(id=tab, seurat=seurat),
                        cluster_picker.ui(id=tab, seurat=seurat, resolution=TRUE, label_switch=TRUE),
                        feature_picker.ui(id=tab, seurat=seurat))
    })  -> server_output$right_sidebar.data_opts
      renderUI({tagList(point_size.ui(id=tab),
                        opacity.ui(id=tab))}) -> server_output$right_sidebar.plotting_opts}})

  # call the modules for this tab
  dimension_reduction <- callModule(module=dimension_reduction.server, id='', seurat=seurat, regex='.*')
  cluster_resolution <- callModule(module=cluster_picker.server, id='', seurat=seurat)
  point_size <- callModule(module=point_size.server, id='')
  opacity <- callModule(module=opacity.server, id='')
  feature_picker <- callModule(module=feature_picker.server, id='', seurat=seurat)
  colour_picker <- list(low='linen', mid='white', high='darkviolet', background=rgb(255, 255, 255, 255, max=255))

  callModule(module=dimension_reduction.show_cluster_idents.server, id='all_clusters', dimension_reduction=dimension_reduction, picked_colours=colour_picker, opacity=opacity, point_size=point_size, cluster_resolution=cluster_resolution)
  callModule(module=dimension_reduction.highlight_feature.server, id='picked_feature', dimension_reduction=dimension_reduction, picked_feature=feature_picker, picked_colours=colour_picker, opacity=opacity, point_size=point_size)
  callModule(module=feature_value_per_cluster.server, id='feature_value_per_cluster', picked_feature=feature_picker, cluster_resolution=cluster_resolution, picked_colours=colour_picker)
}

highlight_feature_and_clusters_tab.server <- function(input, output, session, server_input, server_output, server_session, seurat) {
  # build the sidebar ui
  tab <- 'highlight_feature_and_clusters_tab-'
  observeEvent(eventExpr=server_input$left_sidebar, handlerExpr={
    if(server_input$left_sidebar=='highlight_feature_and_clusters_tab') {
      shinyjs::addClass(selector='body', class='control-sidebar-open')
      renderUI({tagList(dimension_reduction.ui(id=tab, seurat=seurat),
                        cluster_picker.ui(id=tab, seurat=seurat, resolution=TRUE, label_switch=TRUE),
                        feature_picker.ui(id=tab, seurat=seurat))
    })  -> server_output$right_sidebar.data_opts
      renderUI({tagList(point_size.ui(id=tab),
                        opacity.ui(id=tab))}) -> server_output$right_sidebar.plotting_opts}})

  # call the modules for this tab
  dimension_reduction <- callModule(module=dimension_reduction.server, id='', seurat=seurat, regex='.*')
  cluster_resolution <- callModule(module=cluster_picker.server, id='', seurat=seurat)
  point_size <- callModule(module=point_size.server, id='')
  opacity <- callModule(module=opacity.server, id='')
  feature_picker <- callModule(module=feature_picker.server, id='', seurat=seurat)
  colour_picker <- list(low='linen', mid='white', high='darkviolet', background=rgb(255, 255, 255, 255, max=255))

  callModule(module=dimension_reduction.show_cluster_idents.server, id='all_clusters', dimension_reduction=dimension_reduction, picked_colours=colour_picker, opacity=opacity, point_size=point_size, cluster_resolution=cluster_resolution)
  callModule(module=dimension_reduction.highlight_feature.server, id='picked_feature', dimension_reduction=dimension_reduction, picked_feature=feature_picker, picked_colours=colour_picker, opacity=opacity, point_size=point_size)
  callModule(module=feature_value_per_cluster.server, id='feature_value_per_cluster', picked_feature=feature_picker, cluster_resolution=cluster_resolution, picked_colours=colour_picker)
}