
shinyAppUI <- function(...) {
  # tab definitions

  ## cell filtering tab
  cell_filtering.tab <- menuItem(text='Cell filtering', tabName='cell_filtering-tab', icon=icon('filter'))

  ### define layout boxes
  cell_filtering.plot_boxes.defaults <- list(solidHeader=TRUE, width=12, collapsible=TRUE, collapsed=FALSE)
  cell_filtering.boxes <- list(total_expression=append(cell_filtering.plot_boxes.defaults,
                                                       list(title='UMI per cell', footer='Total UMI attributed to a cell', status='primary',
                                                            column(width=2, offset=1, boxplot_plot.ui(id='cell_filtering', feature='nCount_RNA')),
                                                            column(width=3, offset=0, knee_plot.ui(id='cell_filtering', feature='nCount_RNA')),
                                                            column(width=4, offset=1, density_plot.ui(id='cell_filtering', feature='nCount_RNA')))) %>%
                                                do.call(what=box),

                               unique_features=append(cell_filtering.plot_boxes.defaults,
                                                      list(title='Features per cell', footer='Number of distinct features detected in a cell', status='primary',
                                                           column(width=2, offset=1, boxplot_plot.ui(id='cell_filtering', feature='nFeature_RNA')),
                                                           column(width=3, offset=0, knee_plot.ui(id='cell_filtering', feature='nFeature_RNA')),
                                                           column(width=4, offset=1, density_plot.ui(id='cell_filtering', feature='nFeature_RNA')))) %>%
                                               do.call(what=box),

                               percent_mitochondria=append(cell_filtering.plot_boxes.defaults,
                                                           list(title='Mitochondrial UMI', footer='Proportion of mitochondrial genes detected in a cell', status='primary',
                                                                column(width=2, offset=1, boxplot_plot.ui(id='cell_filtering', feature='percent_mt')),
                                                                column(width=3, offset=0, knee_plot.ui(id='cell_filtering', feature='percent_mt')),
                                                                column(width=4, offset=1, density_plot.ui(id='cell_filtering', feature='percent_mt')))) %>%
                                                    do.call(what=box),

                               thresholds=append(cell_filtering.plot_boxes.defaults,
                                                 list(title='Thresholds', status='success',
                                                      column(width=2, total_umi_per_cell_filter.ui(id='cell_filtering')),
                                                      column(width=2, features_per_cell_filter.ui(id='cell_filtering')),
                                                      column(width=2, percent_mt_per_cell_filter.ui(id='cell_filtering')),
                                                      column(width=4, show_filtering_parameters.ui(id='cell_filtering')))) %>%
                                          do.call(what=box))

  ### assemble tab content
  cell_filtering.content <- tabItem(tabName='cell_filtering-tab',
                                    h1('Identify cells that can be removed with quality filters'),
                                    fluidRow(project_name_text_box.ui(id='cell_filtering', width=4),
                                             number_of_reads_text_box.ui(id='cell_filtering', width=2),
                                             number_of_cells_text_box.ui(id='cell_filtering', width=2),
                                             number_of_reads_per_cell_text_box.ui(id='cell_filtering', width=2, f=median),
                                             number_of_genes_per_cell_text_box.ui(id='cell_filtering', width=2, f=median)),
                                    fluidRow(cell_filtering.boxes$total_expression),
                                    fluidRow(cell_filtering.boxes$unique_features),
                                    fluidRow(cell_filtering.boxes$percent_mitochondria),
                                    fluidRow(cell_filtering.boxes$thresholds))

  ## gene highlighting on map tab
  tab <- 'gene_highlighting'
  gene_highlighting.tab <- menuItem(text='Highlight features', tabName=NS(namespace=tab, id='tab'), icon=icon('search'))

  ### define ui elements

  ### define layout boxes
  gene_highlighting.boxes <- list(cluster_dim_plot=box(title='Clustered map',
                                                       footer='Map showing cell clusters',
                                                       status='primary', solidHeader=TRUE, width=4, collapsible=TRUE,
                                                       reduced_dimension_plot.ui(id=tab, feature='selected_cluster_resolution')),
                                  gene_expression_map=box(title='Feature on the map',
                                                          footer='Cells coloured by selected feature',
                                                          status='primary', solidHeader=TRUE, width=4, collapsible=TRUE,
                                                          reduced_dimension_plot.ui(id=tab, feature='picked_feature_values')),
                                  expression_per_cluster=box(title='Feature in clusters',
                                                             footer='Feature values in clusters',
                                                             status='primary', solidHeader=TRUE, width=4, collapsible=TRUE,
                                                             feature_values_per_cluster_plot.ui(id=tab)),
                                  gene_selector=box(title='Select feature',
                                                    status='success', solidHeader=TRUE, width=4, collapsible=TRUE,
                                                    feature_picker.ui(id=tab, include_metadata_switch=TRUE),
                                                    cluster_resolution_picker.ui(id=tab, include_label_switch=TRUE),
                                                    reduction_method_picker.ui(id=tab),
                                                    assay_picker.ui(id=tab)),
                                  plot_options=box(title='Plot options',
                                                   status='success', solidHeader=TRUE, width=4, collapsible=TRUE,
                                                   colour_palette.ui(id=tab, include_full=TRUE,
                                                                     selectors=list(list(inputId='low', label='Low', value='linen'),
                                                                                    list(inputId='high', label='High', value='darkviolet'))),
                                                   point_size_slider.ui(id=tab),
                                                   opacity_slider.ui(id=tab)))

  ### assemble tab content
  gene_highlighting.content <- tabItem(tabName=NS(namespace=tab, id='tab'),
                                       h1('Highlight cell features on the map'),
                                       fluidRow(project_name_text_box.ui(id=tab, width=5),
                                                picked_feature_and_description_text_box.ui(id=tab, width=7)),
                                       fluidRow(number_of_reads_text_box.ui(id=tab, width=2),
                                                number_of_clusters_text_box.ui(id=tab, width=2),
                                                number_of_cells_text_box.ui(id=tab, width=2),
                                                number_of_genes_in_assay_text_box.ui(id=tab, width=2),
                                                number_of_reads_per_cell_text_box.ui(id=tab, width=2, f=median),
                                                number_of_genes_per_cell_text_box.ui(id=tab, width=2, f=median)),
                                       fluidRow(gene_highlighting.boxes$cluster_dim_plot,
                                                gene_highlighting.boxes$gene_expression_map,
                                                gene_highlighting.boxes$expression_per_cluster),
                                       fluidRow(gene_highlighting.boxes$gene_selector,
                                                gene_highlighting.boxes$plot_options))

  ## gene highlighting on map tab
  tab <- 'highlight_features_clusters'
  highlight_features_clusters.tab <- menuItem(text='Highlight features & clusters', tabName=NS(namespace=tab, id='tab'), icon=icon('search'))

  ### define ui elements

  ### define layout boxes
  highlight_features_clusters.boxes <- list(cluster_dim_plot=box(title='Clustered map',
                                                                  footer='Map showing cell clusters',
                                                                  status='primary', solidHeader=TRUE, width=3, collapsible=TRUE,
                                                                  reduced_dimension_plot.ui(id=tab, feature='selected_cluster_resolution')),
                                            sel_cluster_dim_plot=box(title='Clustered map',
                                                                     footer='Map showing selected cell clusters',
                                                                     status='primary', solidHeader=TRUE, width=3, collapsible=TRUE,
                                                                     reduced_dimension_plot.ui(id=tab, feature='selected_cluster_ids')),
                                             gene_expression_map=box(title='Feature on the map',
                                                                     footer='Cells coloured by selected feature',
                                                                     status='primary', solidHeader=TRUE, width=3, collapsible=TRUE,
                                                                     reduced_dimension_plot.ui(id=tab, feature='picked_feature_values')),
                                             expression_per_cluster=box(title='Feature in clusters',
                                                                        footer='Feature values in clusters',
                                                                        status='primary', solidHeader=TRUE, width=3, collapsible=TRUE,
                                                                        feature_values_per_cluster_plot.ui(id=tab)),
                                             gene_selector=box(title='Select feature',
                                                               status='success', solidHeader=TRUE, width=4, collapsible=TRUE,
                                                               feature_picker.ui(id=tab, include_metadata_switch=TRUE),
                                                               cluster_resolution_picker.ui(id=tab, include_label_switch=TRUE),
                                                               cluster_id_picker.ui(id=tab),
                                                               reduction_method_picker.ui(id=tab),
                                                               assay_picker.ui(id=tab)),
                                             plot_options=box(title='Plot options',
                                                              status='success', solidHeader=TRUE, width=4, collapsible=TRUE,
                                                              colour_palette.ui(id=tab, include_full=TRUE,
                                                                                selectors=list(list(inputId='low', label='Low', value='linen'),
                                                                                               list(inputId='high', label='High', value='darkviolet'))),
                                                              point_size_slider.ui(id=tab),
                                                              opacity_slider.ui(id=tab)))

  ### assemble tab content
  highlight_features_clusters.content <- tabItem(tabName=NS(namespace=tab, id='tab'),
                                       h1('Highlight cell features on the map'),
                                       fluidRow(project_name_text_box.ui(id=tab, width=5),
                                                picked_feature_and_description_text_box.ui(id=tab, width=7)),
                                       fluidRow(number_of_reads_text_box.ui(id=tab, width=2),
                                                number_of_clusters_text_box.ui(id=tab, width=2),
                                                number_of_cells_text_box.ui(id=tab, width=2),
                                                number_of_genes_in_assay_text_box.ui(id=tab, width=2),
                                                number_of_reads_per_cell_text_box.ui(id=tab, width=2, f=median),
                                                number_of_genes_per_cell_text_box.ui(id=tab, width=2, f=median)),
                                       fluidRow(highlight_features_clusters.boxes$cluster_dim_plot,
                                                highlight_features_clusters.boxes$sel_cluster_dim_plot,
                                                highlight_features_clusters.boxes$gene_expression_map,
                                                highlight_features_clusters.boxes$expression_per_cluster),
                                       fluidRow(highlight_features_clusters.boxes$gene_selector,
                                                highlight_features_clusters.boxes$plot_options))

  ## submit/configure data tab
  submit_data.tab <- menuItem(text='Configure', tabName='submit_data_tab', icon=icon('cogs'), badgeLabel='!', badgeColor='yellow', selected=TRUE)

  ### assemble tab content
  load_dataset.content <- tabItem(tabName='submit_data_tab',
                                  h1('Select a Seurat object'),
                                  available_seurats.ui(id='load_dataset'))

  ## menu tab hyperlinks
  email_me.tab <- menuItem(text='mail me', href='mailto:christopher.barrington@crick.ac.uk?subject=[seuratvis] Hello there', icon=icon('comment-dots'), newtab=FALSE)
  bug_report.tab <- menuItem(text='report a bug', href='mailto:christopher.barrington@crick.ac.uk?subject=[seuratvis] I found a bug!', icon=icon('bug'), newtab=FALSE)
  feature_request.tab <- menuItem(text='suggest a feature', href='mailto:christopher.barrington@crick.ac.uk?subject=[seuratvis] It would be cool if...', icon=icon('lightbulb'), newtab=FALSE)
  github_link.tab <- menuItem(text=sprintf('GitHub (version: %s)', packageVersion('seuratvis')), href='github.com/ChristopherBarrington/seuratvis', icon=icon('code-branch'), newtab=FALSE)

  # header definition
  logo <- htmltools::HTML("<p style='font-size:26px'>seurat<b>vis</b></p>")
  dashboard_header <- dashboardHeader(title=logo)

  # dashboard body definition
  css <- 'table.dataTable tr.active td, table.dataTable td.active {background-color: #3C8DBC !important;}'
  list(cell_filtering.content,
       gene_highlighting.content,
       highlight_features_clusters.content,
       load_dataset.content) %>%
    do.call(what=tabItems) %>%
    dashboardBody(rclipboardSetup(), tags$head(tags$style(HTML(css)))) -> dashboard_body

  # sidebar definition
  list(cell_filtering.tab,
       gene_highlighting.tab,
       highlight_features_clusters.tab,
       submit_data.tab,
       hr(),
       email_me.tab,
       feature_request.tab,
       bug_report.tab,
       github_link.tab) %>%
    sidebarMenu(id='sidebarmenu') %>%
    dashboardSidebar() -> dashboard_sidebar

  # assemble the final UI
  list(header=dashboard_header,
       sidebar=dashboard_sidebar,
       body=dashboard_body,
       title='seuratvis',
       skin='blue') %>%
    do.call(what=dashboardPage)
}
