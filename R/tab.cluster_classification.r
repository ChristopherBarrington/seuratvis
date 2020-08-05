#'
#' 
cluster_classification.tab <- function() {
  bquote({
    menuItem(text='Cluster classification', icon=icon('object-group'), startExpanded=TRUE,
             menuSubItem(text='FindMarkers results', tabName='findmarkers_results_tab', icon=menuSubItem_icon()),
             menuItem(text='Gene modules', startExpanded=TRUE, icon=icon('layer-group'),
                      menuSubItem(text='Score in clusters', tabName='gene_module_score_in_clusters_tab', icon=menuSubItem_icon()),
                      menuSubItem(text='Scores in a cluster', tabName='gene_module_scores_in_a_cluster_tab', icon=menuSubItem_icon())) %>%
              modify_stop_propagation()) %>%
      modify_stop_propagation() -> menu_item

    list(tabItem(tabName='findmarkers_results_tab',
                 h1('Feature markers of clusters'),
                 fluidRow(dataset_info_text_box.ui(id=NS('findmarkers_results_tab', 'project_name'), width=12)),
                 fluidRow(boxPlus(title='Results table', closable=FALSE, width=12, status='primary'), DT::dataTableOutput(outputId=NS('findmarkers_results_tab', 'table')) %>% withSpinner())),
         
         tabItem(tabName='gene_module_score_in_clusters_tab',
                 h1('Gene module score in clusters'),
                 fluidRow(dataset_info_text_box.ui(id=NS('gene_module_score_in_clusters_tab', 'project_name'), width=12)),
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
                 fluidRow(dataset_info_text_box.ui(id=NS('gene_module_scores_in_a_cluster_tab', 'project_name'), width=12)),

                 # fluidRow(column(width=12, feature_ridges.plot(id=NS('gene_module_scores_in_a_cluster_tab', 'scores_plot'))))

                 column(width=7,
                        boxPlus(title='Gene module score', closable=FALSE, width=12, height='75vh', status='primary', 
                                feature_ridges.plot(id=NS('gene_module_scores_in_a_cluster_tab', 'scores_plot')))),
                 column(width=5,
                        boxPlus(title='Map', closable=FALSE, width=12, height='35vh', status='primary',
                                dimension_reduction.plot(id=NS('gene_module_scores_in_a_cluster_tab', 'map'))),
                        boxPlus(title='Genes', closable=FALSE, width=12, height='35vh', status='primary',
                                genes_in_modules.table(id=NS('gene_module_scores_in_a_cluster_tab', 'modules_table'))))

                 )) -> content

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
      if(nrow(seurat$FindMarkersResults$table)==0) {
        error_alert(session=session, title='FindMarkers results', text='This Seurat object does not have any FindMarkers results.')
        go_to_config(session=server_session)
        return(NULL)
      }

      # tab %<>% str_c('-')
      renderUI({tagList(filterDF_UI(id=NS(tab, 'filter_parameters')))})  -> server_output$right_sidebar.data_opts
      renderUI({tagList()}) -> server_output$right_sidebar.plotting_opts

      sendSweetAlert(session=session, type='warning',
                     title='Notice', btn_labels='OK',
                     text='For some reason you have to reload the Seurat object from the config tab in the right sidebar before the table will render. This will be fixed later.',
                     closeOnClickOutside=TRUE, showCloseButton=FALSE)
    }})

  # call the modules for this tab
  callModule(module=dataset_info_text_box.project_name, id='project_name', seurat=seurat)

# if(!is.null(reactive(seurat$FindMarkersResults))) {
  # handle the data.frame filtering
  ## call the data.frame filtering filtering module
  callModule(session=session, module=filterDF, id='filter_parameters',
             picker=TRUE, drop_ids=FALSE,
             data_name=reactive(seurat$formatted_project_name), 
             data_table=reactive(seurat$FindMarkersResults$table),
             data_vars=reactive(seurat$FindMarkersResults$vars)) -> results_filter

  DT::renderDataTable({
    DT::datatable(data=results_filter$data_filtered(),
                  rownames=FALSE,
                  options=list(columnDefs=list(list(className='dt-right', targets=c(1,6))),
                               ordering=FALSE,
                               dom='litp'),
                  style='bootstrap4',
                  class='stripe') %>%
    formatRound(columns=c('Cluster detection', 'Map detection'), digits=2) %>%
    formatRound(columns=c('Avg. logFC'), digits=3)}) -> output$table
# }
}

#'
#' 
gene_module_score_in_clusters_tab.server <- function(input, output, session, server_input, server_output, server_session, seurat) {
  # build the sidebar ui
  observeEvent(eventExpr=server_input$left_sidebar, handlerExpr={
    tab <- 'gene_module_score_in_clusters_tab'
    if(server_input$left_sidebar==tab) {
      tab %<>% str_c('-')
      renderUI({tagList(cluster_picker.ui(id=tab, seurat=seurat, resolution=TRUE, picker=TRUE, label_switch=FALSE),
                        feature_picker.ui(id=tab, seurat=seurat, choices=list(`Gene modules`='gene_modules'), selected='gene_modules', gene_modules_opts=list(multiple=FALSE), include_feature_type=FALSE, include_values_range=FALSE),
                        dimension_reduction.ui(id=tab, seurat=seurat))})  -> server_output$right_sidebar.data_opts
      renderUI({tagList()}) -> server_output$right_sidebar.plotting_opts}})

  # call the modules for this tab
  dimension_reduction <- callModule(module=dimension_reduction.server, id='', seurat=seurat)
  cluster_resolution <- callModule(module=cluster_picker.server, id='', seurat=seurat)
  feature_picker <- callModule(module=feature_picker.server, id='', seurat=seurat)
  colour_picker <- list(low='linen', mid='white', high='darkviolet', background=rgb(255, 255, 255, 255, max=255))

  callModule(module=dataset_info_text_box.project_name, id='project_name', seurat=seurat)
  callModule(module=feature_ridge_by_idents.server, id='scores_plot', picked_feature=feature_picker, picked_clusters=cluster_resolution)
  callModule(module=dimension_reduction.show_selected_clusters.server, id='map', dimension_reduction=dimension_reduction, point_size=list(size=0.6), cluster_resolution=cluster_resolution, picked_colours=colour_picker)
  callModule(module=genes_in_modules.server, id='modules_table', seurat=seurat, picked_feature=feature_picker)
}

#'
#' 
gene_module_scores_in_a_cluster_tab.server <- function(input, output, session, server_input, server_output, server_session, seurat) {
  # build the sidebar ui
  observeEvent(eventExpr=server_input$left_sidebar, handlerExpr={
    tab <- 'gene_module_scores_in_a_cluster_tab'
    if(server_input$left_sidebar==tab) {
      tab %<>% str_c('-')
      renderUI({tagList(cluster_picker.ui(id=tab, seurat=seurat, resolution=TRUE, picker=TRUE, label_switch=FALSE, multi_picker=FALSE),
                        feature_picker.ui(id=tab, seurat=seurat, choices=list(`Gene modules`='gene_modules'), selected='gene_modules', gene_modules_opts=list(multiple=TRUE), include_feature_type=FALSE, include_values_range=FALSE),
                        dimension_reduction.ui(id=tab, seurat=seurat))})  -> server_output$right_sidebar.data_opts
      renderUI({tagList()}) -> server_output$right_sidebar.plotting_opts}})

  # call the modules for this tab
  dimension_reduction <- callModule(module=dimension_reduction.server, id='', seurat=seurat)
  cluster_resolution <- callModule(module=cluster_picker.server, id='', seurat=seurat)
  feature_picker <- callModule(module=feature_picker.server, id='', seurat=seurat)
  colour_picker <- list(low='linen', mid='white', high='darkviolet', background=rgb(255, 255, 255, 255, max=255))

  callModule(module=dataset_info_text_box.project_name, id='project_name', seurat=seurat)
  callModule(module=feature_ridges_by_ident.server, id='scores_plot', picked_feature=feature_picker, picked_clusters=cluster_resolution)
  callModule(module=dimension_reduction.show_selected_clusters.server, id='map', dimension_reduction=dimension_reduction, point_size=list(size=0.6), cluster_resolution=cluster_resolution, picked_colours=colour_picker)
  callModule(module=genes_in_modules.server, id='modules_table', seurat=seurat, picked_feature=feature_picker)
}

