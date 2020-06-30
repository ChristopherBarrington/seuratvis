#'
#' 
preprocessing.tab <- function() {
  bquote({
    # define the menu for this group of tabs
    menuItem(text='Preprocessing', icon=icon('toolbox'), startExpanded=TRUE,
             menuSubItem(text='Cell filtering', tabName='cell_filtering_tab', icon=menuSubItem_icon()),
             menuSubItem(text='Dimensionality', tabName='dimensionality_tab', icon=menuSubItem_icon()),
             menuSubItem(text='Cluster filtering', tabName='cluster_filtering_tab', icon=under_construction_icon())) %>%
      modify_stop_propagation() -> menu_item

    # define the content of each tab
    ## filter cells by per-cell metrics
    tab_name <- 'cell_filtering_tab'
    tabItem(tabName=tab_name, h1('Identify cells that can be removed with quality filters'),
            fluidRow(dataset_info_text_box.ui(id=NS(tab_name, 'project_name'), width=12),
                     dataset_info_text_box.ui(id=NS(tab_name, 'n_umi'), width=3),
                     dataset_info_text_box.ui(id=NS(tab_name, 'n_cells'), width=3),
                     dataset_info_text_box.ui(id=NS(tab_name, 'n_umi_per_cell'), width=3),
                     dataset_info_text_box.ui(id=NS(tab_name, 'n_features_per_cell'), width=3)),
            fluidRow(boxPlus(title='UMI per cell', closable=FALSE, width=12, status='primary',
                             column(width=2, offset=1, boxplot.plot(id=NS(tab_name, 'n_umi'))),
                             column(width=3, offset=0, knee.plot(id=NS(tab_name, 'n_umi'))),
                             column(width=4, offset=1, density.plot(id=NS(tab_name, 'n_umi'))))),
            fluidRow(boxPlus(title='Features per cell', closable=FALSE, width=12, status='primary',
                             column(width=2, offset=1, boxplot.plot(id=NS(tab_name, 'n_features'))),
                             column(width=3, offset=0, knee.plot(id=NS(tab_name, 'n_features'))),
                             column(width=4, offset=1, density.plot(id=NS(tab_name, 'n_features'))))),
            fluidRow(boxPlus(title='Mitochondrial UMI', closable=FALSE, width=12, status='primary',
                             column(width=2, offset=1, boxplot.plot(id=NS(tab_name, 'proportion_mt'))),
                             column(width=3, offset=0, knee.plot(id=NS(tab_name, 'proportion_mt'))),
                             column(width=4, offset=1, density.plot(id=NS(tab_name, 'proportion_mt')))))) -> cell_filtering_tab
    ## dimensionality tab
    tab_name <- 'dimensionality_tab'
    tabItem(tabName='dimensionality_tab', h1('Determine the dimensionality of a dataset'),
            fluidRow(dataset_info_text_box.ui(id=NS('dimensionality_tab', 'project_name'), width=12)),
            fluidRow(boxPlus(title='Elbow plot', closable=FALSE, width=6, status='primary', dimensionality.plot(id=NS('dimensionality_tab','elbow'))),
                     boxPlus(title='Top contributing features', closable=FALSE, width=6, status='primary', dimensionality.plot(id=NS('dimensionality_tab','top_features'))),
                     boxPlus(title='JackStraw plot', closable=FALSE, width=6, status='primary', dimensionality.plot(id=NS('dimensionality_tab','jackstraw'))),
                     boxPlus(title='JackStraw p-values', closable=FALSE, width=6, status='primary', dimensionality.plot(id=NS('dimensionality_tab','jackstraw_pvalue'))))) -> dimensionality_tab
    
    ## filter clusters using cluster-wide metrics
    tab_name <- ns <- 'cluster_filtering_tab'
    ns %<>% NS()
    tabItem(tabName=tab_name, h1('Remove clusters of cells before reanalysis'),
            p('Use these tools to identify clusters from the reduced dimension map that could be removed. Clusters that represent contaminant or are characterised by high proportion of mitochondrial RNA, for example, may be suitable for removal before the remaining cell set is re-analysed.'),
            fluidRow(dataset_info_text_box.ui(id=ns('project_name'), width=7),
                     picked_feature_and_description_text_box.ui(id=ns('feature_description'), width=5),
                     dataset_info_text_box.ui(id=ns('n_umi'), width=3),
                     dataset_info_text_box.ui(id=ns('n_cells'), width=3),
                     dataset_info_text_box.ui(id=ns('n_umi_per_cell'), width=3),
                     dataset_info_text_box.ui(id=ns('n_features_per_cell'), width=3)),
            fluidRow(boxPlus(title='Selected feature', closable=FALSE, width=4, dimension_reduction.plot(id=ns('picked_feature'))),
                     boxPlus(title='Selected cluster(s)', closable=FALSE, width=4, dimension_reduction.plot(id=ns('picked_clusters'))),
                     boxPlus(title='Feature values in clusters', closable=FALSE, width=4, feature_value_per_cluster.plot(id=ns('feature_value_per_cluster'))))) -> cluster_filtering_tab

    # collect all of the contents for this tab group
    content <- list(cell_filtering_tab, dimensionality_tab, cluster_filtering_tab)

    # add this menu and content to the app
    menus %<>% append(list(menu_item))
    contents %<>% append(content)})
}

#'
#' 
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
  n_umi_output <- callModule(module=density_plot.n_umi.server, id='n_umi', seurat=seurat)
  n_features_output <- callModule(module=density_plot.n_features.server, id='n_features', seurat=seurat)
  proportion_mt_output <- callModule(module=density_plot.proportion_mt.server, id='proportion_mt', seurat=seurat)

  filtering_parameters <- reactiveValues()
  n_umi_filter <- callModule(module=filter_n_umi.server, id='', seurat=seurat, cell_filtering=filtering_parameters, linked_density_plot=n_umi_output)
  n_features_filter <- callModule(module=filter_n_features.server, id='', seurat=seurat, cell_filtering=filtering_parameters, linked_density_plot=n_features_output)
  proportion_mt_filter <- callModule(module=filter_proportion_mt.server, id='', seurat=seurat, cell_filtering=filtering_parameters, linked_density_plot=proportion_mt_output)
  filtering_parameters <- callModule(module=dataset_filtering.server, id='filtering', seurat=seurat, filters=list(n_umi_filter, n_features_filter, proportion_mt_filter))

  callModule(module=dataset_info_text_box.project_name, id='project_name', seurat=seurat)
  callModule(module=dataset_info_text_box.n_filtered_umi, id='n_umi', seurat=seurat, cell_filtering=filtering_parameters)
  callModule(module=dataset_info_text_box.n_filtered_cells, id='n_cells', seurat=seurat, cell_filtering=filtering_parameters)
  callModule(module=dataset_info_text_box.n_umi_per_filtered_cell, id='n_umi_per_cell', seurat=seurat, cell_filtering=filtering_parameters)
  callModule(module=dataset_info_text_box.n_features_per_filtered_cell, id='n_features_per_cell', seurat=seurat, cell_filtering=filtering_parameters)
  callModule(module=show_filtering_parameters.server, id='', seurat=seurat, cell_filtering=filtering_parameters, filters=list(n_umi_filter, n_features_filter, proportion_mt_filter))

  callModule(module=boxplot_plot.n_features.server, id='n_features', seurat=seurat, cell_filtering=filtering_parameters)
  callModule(module=boxplot_plot.n_umi.server, id='n_umi', seurat=seurat, cell_filtering=filtering_parameters)
  callModule(module=boxplot_plot.proportion_mt.server, id='proportion_mt', seurat=seurat, cell_filtering=filtering_parameters)

  callModule(module=knee_plot.n_features.server, id='n_features', seurat=seurat, cell_filtering=filtering_parameters)
  callModule(module=knee_plot.n_umi.server, id='n_umi', seurat=seurat, cell_filtering=filtering_parameters)
  callModule(module=knee_plot.proportion_mt.server, id='proportion_mt', seurat=seurat, cell_filtering=filtering_parameters)

}

#'
#' 
dimensionality_tab.server <- function(input, output, session, server_input, server_output, server_session, seurat) {
  # build the sidebar ui
  observeEvent(eventExpr=server_input$left_sidebar, handlerExpr={
    tab <- 'dimensionality_tab'
    if(server_input$left_sidebar==tab) {
      tab %<>% str_c('-')
      renderUI({tagList(dimension_reduction.ui(id=tab, seurat=seurat, regex='^pca'),
                        component_picker_slider.ui(id=tab))})  -> server_output$right_sidebar.data_opts
      renderUI({tagList()}) -> server_output$right_sidebar.plotting_opts}})

  # call the modules for this tab
  reduction_picker <- callModule(module=dimension_reduction.server, id='', seurat=seurat)
  components_picker <- callModule(module=components_selector.server, id='', seurat=seurat, picked_reduction=reduction_picker)

  callModule(module=dataset_info_text_box.project_name, id='project_name', seurat=seurat)

  callModule(module=dimensionality.elbow, id='elbow', seurat=seurat, picked_reduction=reduction_picker)
  callModule(module=dimensionality.jackstraw, id='jackstraw', seurat=seurat, picked_reduction=reduction_picker, picked_components=components_picker)
  callModule(module=dimensionality.jackstraw_pvalue, id='jackstraw_pvalue', seurat=seurat, picked_reduction=reduction_picker)
  callModule(module=dimensionality.top_features_pca_heatmap, id='top_features', seurat=seurat, picked_reduction=reduction_picker, picked_components=components_picker)
}

#'
#' 
cluster_filtering_tab.server <- function(input, output, session, server_input, server_output, server_session, seurat) {
  tab <- session$ns('') %>% str_remove('-$')
  
  # build the sidebar ui
  observeEvent(eventExpr=server_input$left_sidebar, handlerExpr={
    if(server_input$left_sidebar==tab) {    
      tab %<>% str_c('-')
      renderUI({tagList(dimension_reduction.ui(id=tab, seurat=seurat),
                        cluster_picker.ui(id=tab, seurat=seurat, resolution=TRUE, picker=FALSE, label_switch=TRUE),
                        feature_picker.ui(id=tab, seurat=seurat),
                        show_filtering_parameters.ui(id=tab, label='Cell filtering parameters'))})  -> server_output$right_sidebar.data_opts
      renderUI({tagList()}) -> server_output$right_sidebar.plotting_opts}})

  # call the modules for this tab
  # n_umi_output <- callModule(module=density_plot.n_umi.server, id='n_umi', seurat=seurat)
  # n_features_output <- callModule(module=density_plot.n_features.server, id='n_features', seurat=seurat)
  # proportion_mt_output <- callModule(module=density_plot.proportion_mt.server, id='proportion_mt', seurat=seurat)

  # filtering_parameters <- reactiveValues()
  # n_umi_filter <- callModule(module=filter_n_umi.server, id='', seurat=seurat, cell_filtering=filtering_parameters, linked_density_plot=n_umi_output)
  # n_features_filter <- callModule(module=filter_n_features.server, id='', seurat=seurat, cell_filtering=filtering_parameters, linked_density_plot=n_features_output)
  # proportion_mt_filter <- callModule(module=filter_proportion_mt.server, id='', seurat=seurat, cell_filtering=filtering_parameters, linked_density_plot=proportion_mt_output)
  # filtering_parameters <- callModule(module=dataset_filtering.server, id='filtering', seurat=seurat, filters=list(n_umi_filter, n_features_filter, proportion_mt_filter))

  dimension_reduction <- callModule(module=dimension_reduction.server, id='', seurat=seurat, regex='.*')
  cluster_resolution <- callModule(module=cluster_picker.server, id='', seurat=seurat)
  feature_picker <- callModule(module=feature_picker.server, id='', seurat=seurat)
  # point_size <- callModule(module=point_size.server, id='')
  # opacity <- callModule(module=opacity.server, id='')
  colour_picker <- list(low='linen', mid='white', high='darkviolet', background=rgb(255, 255, 255, 255, max=255))
  point_size <- list(size=0.6)
  opacity <- list(alpha=1)

  callModule(module=dataset_info_text_box.project_name, id='project_name', seurat=seurat)
  callModule(module=picked_feature_and_description_text_box.server, id='feature_description', seurat=seurat, picked_feature=feature_picker)
  callModule(module=dataset_info_text_box.n_filtered_umi, id='n_umi', seurat=seurat, cell_filtering=filtering_parameters)
  callModule(module=dataset_info_text_box.n_filtered_cells, id='n_cells', seurat=seurat, cell_filtering=filtering_parameters)
  callModule(module=dataset_info_text_box.n_umi_per_filtered_cell, id='n_umi_per_cell', seurat=seurat, cell_filtering=filtering_parameters)
  callModule(module=dataset_info_text_box.n_features_per_filtered_cell, id='n_features_per_cell', seurat=seurat, cell_filtering=filtering_parameters)

  callModule(module=show_filtering_parameters.server, id='', seurat=seurat, cell_filtering=filtering_parameters)

  callModule(module=dimension_reduction.highlight_feature.server, id='picked_feature', dimension_reduction=dimension_reduction, picked_feature=feature_picker, picked_colours=colour_picker, opacity=opacity, point_size=point_size)
  callModule(module=dimension_reduction.show_selected_clusters.server, id='picked_clusters', dimension_reduction=dimension_reduction, picked_colours=colour_picker, point_size=point_size, cluster_resolution=cluster_resolution)
  callModule(module=feature_value_per_cluster.server, id='feature_value_per_cluster', picked_feature=feature_picker, cluster_resolution=cluster_resolution, picked_colours=colour_picker)

}
