
highlight_features.tab <- function() {
  bquote({
    menuItem(text='Highlight features', icon=icon('highlighter'), startExpanded=TRUE,
             menuSubItem(text='Highlight feature', tabName='highlight_feature_tab', icon=menuSubItem_icon()),
             menuSubItem(text='Highlight feature and cluster', tabName='highlight_feature_and_clusters_tab', icon=menuSubItem_icon()),
             menuSubItem(text='Highlight multiple features', tabName='highlight_multiple_features_tab', icon=icon('tools'))) %>%
    modify_stop_propagation() -> menu_item

    list(tabItem(tabName='highlight_feature_tab',
                 h1('Highlight a feature on a map'),
                 fluidRow(dataset_info_text_box.ui(id=NS('highlight_feature_tab', 'project_name'), width=6),
                          picked_feature_and_description_text_box.ui(id=NS('highlight_feature_tab', 'feature_description'), width=6)),
                 fluidRow(boxPlus(title='All clusters in map', closable=FALSE, width=4, dimension_reduction.plot(id='highlight_feature_tab-all_clusters')),
                          boxPlus(title='Selected feature', closable=FALSE, width=4, dimension_reduction.plot(id='highlight_feature_tab-picked_feature')),
                          boxPlus(title='Feature values in clusters', closable=FALSE, width=4, feature_value_per_cluster.plot(id='highlight_feature_tab-feature_value_per_cluster')))),
         
         tabItem(tabName='highlight_feature_and_clusters_tab',
                 h1('Highlight a feature and clusters'),
                 fluidRow(dataset_info_text_box.ui(id=NS('highlight_feature_and_clusters_tab', 'project_name'), width=6),
                          picked_feature_and_description_text_box.ui(id=NS('highlight_feature_and_clusters_tab', 'feature_description'), width=6)),
                 fluidRow(boxPlus(title='All clusters in map', closable=FALSE, width=3, dimension_reduction.plot(id='highlight_feature_and_clusters_tab-all_clusters')),
                          boxPlus(title='Selected feature', closable=FALSE, width=3, dimension_reduction.plot(id='highlight_feature_and_clusters_tab-picked_feature')),
                          boxPlus(title='Selected cluster(s)', closable=FALSE, width=3, dimension_reduction.plot(id='highlight_feature_and_clusters_tab-picked_clusters')),
                          boxPlus(title='Feature values in clusters', closable=FALSE, width=3, feature_value_per_cluster.plot(id='highlight_feature_and_clusters_tab-feature_value_per_cluster')))),

         tabItem(tabName='highlight_multiple_features_tab',
                 h1('Add multiple maps and plot features in cells'),
                 dimension_reduction.plotbox(id='highlight_multiple_features_tab-feature0'))) -> content

    menus %<>% append(list(menu_item))
    contents %<>% append(content)})
}

highlight_feature_tab.server <- function(input, output, session, server_input, server_output, server_session, seurat) {
  # build the sidebar ui
  observeEvent(eventExpr=server_input$left_sidebar, handlerExpr={
    tab <- 'highlight_feature_tab'
    if(server_input$left_sidebar==tab) {    
      tab %<>% str_c('-')
      renderUI({tagList(dimension_reduction.ui(id=tab, seurat=seurat),
                        cluster_picker.ui(id=tab, seurat=seurat, resolution=TRUE, picker=FALSE, label_switch=TRUE),
                        feature_picker.ui(id=tab, seurat=seurat))})  -> server_output$right_sidebar.data_opts
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

  callModule(module=dataset_info_text_box.project_name, id='project_name', seurat=seurat)
  callModule(module=picked_feature_and_description_text_box.server, id='feature_description', seurat=seurat, picked_feature=feature_picker)
}

highlight_feature_and_clusters_tab.server <- function(input, output, session, server_input, server_output, server_session, seurat) {
  # build the sidebar ui
  observeEvent(eventExpr=server_input$left_sidebar, handlerExpr={
    tab <- 'highlight_feature_and_clusters_tab'
    if(server_input$left_sidebar==tab) {
      tab %<>% str_c('-')
      renderUI({tagList(dimension_reduction.ui(id=tab, seurat=seurat),
                        cluster_picker.ui(id=tab, seurat=seurat, resolution=TRUE, picker=TRUE, label_switch=TRUE),
                        feature_picker.ui(id=tab, seurat=seurat))})  -> server_output$right_sidebar.data_opts
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
  callModule(module=dimension_reduction.show_selected_clusters.server, id='picked_clusters', dimension_reduction=dimension_reduction, picked_colours=colour_picker, point_size=point_size, cluster_resolution=cluster_resolution)
  callModule(module=feature_value_per_cluster.server, id='feature_value_per_cluster', picked_feature=feature_picker, cluster_resolution=cluster_resolution, picked_colours=colour_picker)

  callModule(module=dataset_info_text_box.project_name, id='project_name', seurat=seurat)
  callModule(module=picked_feature_and_description_text_box.server, id='feature_description', seurat=seurat, picked_feature=feature_picker)
}

highlight_multiple_features.server <- function(input, output, session, server_input, server_output, server_session, seurat) {
  # build the sidebar ui
  observeEvent(eventExpr=server_input$left_sidebar, handlerExpr={
    tab <- 'highlight_multiple_features_tab'
    if(server_input$left_sidebar==tab) {    
      tab %<>% str_c('-')
      renderUI({tagList(dimension_reduction.ui(id=tab, seurat=seurat),
                        actionBttn(inputId=NS('highlight_multiple_features_tab','add_plot_ui'), label='Add a map',
                                   style='bordered', color='primary', icon=icon('plus'))
                        )})  -> server_output$right_sidebar.data_opts
      renderUI({tagList(point_size.ui(id=tab),
                        opacity.ui(id=tab))}) -> server_output$right_sidebar.plotting_opts

      renderUI({tagList(feature_picker.ui(id=str_c(tab, 'feature0'), seurat=seurat))}) -> output[['feature0-dropdown']]
    }})

  # react to the action button being pressed
  observeEvent(input$add_plot_ui, {
    feature_n <- input$add_plot_ui
    feature_id <- str_c('feature', feature_n)
    box_id <- NS('highlight_multiple_features_tab', feature_id)
    div_id <- sprintf(fmt='#boxid%s', feature_n-1)
    ddn_id <- NS(feature_id, 'dropdown')

    dimension_reduction.plotbox(id=box_id, n=feature_n) -> ui
    insertUI(selector=div_id, where='afterEnd', ui=ui)
    renderUI({tagList(feature_picker.ui(id=box_id, seurat=seurat))}) -> output[[ddn_id]]
 
    feature_pickers[[feature_id]] <- callModule(module=feature_picker.server, id=feature_id, seurat=seurat)
    callModule(module=dimension_reduction.highlight_feature.server, id=feature_id, dimension_reduction=dimension_reduction, picked_feature=feature_pickers[[feature_id]], picked_colours=colour_picker, opacity=opacity, point_size=point_size)

  })

  # call the modules for this tab
  opacity <- callModule(module=opacity.server, id='')
  point_size <- callModule(module=point_size.server, id='')
  colour_picker <- list(low='linen', mid='white', high='darkviolet', background=rgb(255, 255, 255, 255, max=255))
  dimension_reduction <- callModule(module=dimension_reduction.server, id='', seurat=seurat)

  feature_pickers <- list(feature0=callModule(module=feature_picker.server, id='feature0', seurat=seurat))
  callModule(module=dimension_reduction.highlight_feature.server, id='feature0', dimension_reduction=dimension_reduction, picked_feature=feature_pickers$feature0, picked_colours=colour_picker, opacity=opacity, point_size=point_size)

}
