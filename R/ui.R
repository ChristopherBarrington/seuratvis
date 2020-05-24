
shinyAppUI <- function(...) {
  # tab definitions

  ## cell filtering tab
  cell_filtering.tab <- menuItem(text='Cell filtering', tabName='cell_filtering-tab', icon=icon('filter'))

  ### define ui elements

  ### define layout boxes
  cell_filtering.plot_boxes.defaults <- list(solidHeader=TRUE, width=12, collapsible=TRUE, collapsed=TRUE)
  cell_filtering.boxes <- list(total_expression=append(cell_filtering.plot_boxes.defaults,
                                                       list(title='Expression per cell', footer='Total number of reads attributed to a cell', status='success',
                                                            column(width=2, offset=1, plotOutput(outputId='cell_filtering-total_expression_boxplot') %>% withSpinner()),
                                                            column(width=3, offset=0, plotOutput(outputId='cell_filtering-total_expression_knee') %>% withSpinner()),
                                                            column(width=4, offset=1, plotOutput(outputId='cell_filtering-total_expression_density', brush=brushOpts(id='total_expression_density.brush', direction='x')) %>% withSpinner()))) %>%
                                                do.call(what=box),

                               unique_features=append(cell_filtering.plot_boxes.defaults,
                                                      list(title='Genes per cell', footer='Number of distinct genes detected in a cell', status='success',
                                                           column(width=2, offset=1, plotOutput(outputId='cell_filtering-unique_genes_boxplot') %>% withSpinner()),
                                                           column(width=3, offset=0, plotOutput(outputId='cell_filtering-unique_genes_knee') %>% withSpinner()),
                                                           column(width=4, offset=1, plotOutput(outputId='cell_filtering-unique_genes_density', brush=brushOpts(id='unique_genes_density.brush', direction='x')) %>% withSpinner()))) %>%
                                               do.call(what=box),

                               percent_mitochondria=append(cell_filtering.plot_boxes.defaults,
                                                           list(title='Mitochondrial expression', footer='Proportion of mitochondrial genes detected in a cell', status='success',
                                                                column(width=2, offset=1, plotOutput(outputId='cell_filtering-percent_mitochondria_boxplot') %>% withSpinner()),
                                                                column(width=4, offset=0, plotOutput(outputId='cell_filtering-percent_mitochondria_knee') %>% withSpinner()),
                                                                column(width=4, offset=1, plotOutput(outputId='cell_filtering-percent_mitochondria_density', brush=brushOpts(id='percent_mitochondria_density.brush', direction='x')) %>% withSpinner()))) %>%
                                                    do.call(what=box),

                               thresholds=append(cell_filtering.plot_boxes.defaults,
                                                 list(title='Thresholds', status='primary',
                                                      column(width=2, total_umi_per_cell_filter.ui('cell_filtering')),
                                                      column(width=2, features_per_cell_filter.ui('cell_filtering')),
                                                      column(width=2, percent_mt_per_cell_filter.ui(id='cell_filtering')),
                                                      column(width=4, show_filtering_parameters.ui('cell_filtering')))) %>%
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
  gene_highlighting.tab <- menuItem(text='Highlight features', tabName='gene_highlighting-tab', icon=icon('search'))

  ### define ui elements

  ### define layout boxes
  gene_highlighting.boxes <- list(cluster_dim_plot=box(title='Clustered map',
                                                       footer='Map showing cell clusters',
                                                       status='success',
                                                       solidHeader=TRUE,
                                                       width=4,
                                                       collapsible=TRUE,
                                                       {plotOutput('genes_highlighting-clustered_map') %>% withSpinner()}),
                                  gene_expression_map=box(title='Feature on map',
                                                          footer='Cells coloured by selected feature',
                                                          status='success',
                                                          solidHeader=TRUE,
                                                          width=4,
                                                          collapsible=TRUE,
                                                          {plotOutput('genes_highlighting-gene_expression_map') %>% withSpinner()}),
                                  expression_per_cluster=box(title='Feature in clusters',
                                                             footer='Feature values in clusters',
                                                             status='success',
                                                             solidHeader=TRUE,
                                                             width=4,
                                                             collapsible=TRUE,
                                                             {plotOutput('genes_highlighting-expression_per_cluster') %>% withSpinner()}),
                                                             # {plotlyOutput('genes_highlighting-expression_per_cluster') %>% withSpinner()}),
                                  gene_selector=box(title='Select feature',
                                                    status='primary',
                                                    solidHeader=TRUE,
                                                    width=4,
                                                    collapsible=TRUE,
                                                    feature_picker.ui(id='gene_highlighting'),
                                                    cluster_resolution_picker.ui(id='gene_highlighting', include_label_switch=TRUE),
                                                    reduction_method_picker.ui(id='gene_highlighting'),
                                                    assay_picker.ui(id='gene_highlighting')),
                                  plot_options=box(title='Plot options',
                                                   status='primary',
                                                   solidHeader=TRUE,
                                                   width=4,
                                                   collapsible=TRUE,
                                                   colour_palette.ui(id='gene_highlighting', include_full=TRUE,
                                                                     selectors=list(list(inputId='low', label='Low', value='linen'),
                                                                                    list(inputId='high', label='High', value='darkviolet'))),
                                                   point_size_slider.ui(id='gene_highlighting'),
                                                   opacity_slider.ui(id='gene_highlighting')))

  ### assemble tab content
  gene_highlighting.content <- tabItem(tabName='gene_highlighting-tab',
                                       h1('Highlight cell features on the map'),
                                       fluidRow(project_name_text_box.ui(id='gene_highlighting', width=5),
                                                gene_name_and_description_text_box.ui(id='gene_highlighting', width=7)),
                                       fluidRow(number_of_reads_text_box.ui(id='gene_highlighting', width=2),
                                                number_of_clusters_text_box.ui(id='gene_highlighting', width=2),
                                                number_of_cells_text_box.ui(id='gene_highlighting', width=2),
                                                number_of_genes_in_assay_text_box.ui(id='gene_highlighting', width=2),
                                                number_of_reads_per_cell_text_box.ui(id='gene_highlighting', width=2, f=median),
                                                number_of_genes_per_cell_text_box.ui(id='gene_highlighting', width=2, f=median)),
                                       fluidRow(gene_highlighting.boxes$cluster_dim_plot,
                                                gene_highlighting.boxes$gene_expression_map,
                                                gene_highlighting.boxes$expression_per_cluster),
                                       fluidRow(gene_highlighting.boxes$gene_selector,
                                                gene_highlighting.boxes$plot_options))

  ## discover markers tab
  discover_markers.tab <- menuItem(text='Discover markers', tabName='discover_markers-tab', icon=icon('binoculars'), badgeLabel='!', badgeColor='red')

  ### define ui elements

  ### define layout boxes

  ### assemble tab content
  discover_markers.content <- tabItem(tabName = 'discover_markers-tab',
                                      h1('Find gene markers for clusters'))

  ## features of interest heatmap tab
  features_heatmap.tab <- menuItem(text='Features heatmap', tabName='features_heatmap-tab', icon=icon('crosshairs'), badgeLabel='!', badgeColor='red')

  ### define ui elements
  features_heatmap.features_of_interest.tb <- textAreaInput(inputId='features_of_interest.tb', label='Feature names', value='', width='100%', height='100%', placeholder='A list of gene names separated by whitespace, semicolon or comma')
  features_heatmap.use_average_value.checkbox <- checkboxInput(inputId='features_heatmap.use_average_value.checkbox', label='Use mean of cells in cluster', value=TRUE)
  features_heatmap.show_feature_names.checkbox <- checkboxInput(inputId='features_heatmap.show_feature_names.checkbox', label='Show feature names', value=FALSE)
  features_heatmap.seurat_cluster_set.dropdown <- selectInput(inputId = 'features_heatmap.seurat_cluster_set.dd', label = 'Seurat cluster definitions', choices = 'seurat_clusters', selected = 'seurat_clusters', multiple = FALSE)

  ### define layout boxes
  features_heatmap.boxes <- list(heatmap_plot=box(title='Heatmap of selected features',
                                                  footer='Normalised expression of genes in clusters',
                                                  status='success',
                                                  solidHeader=TRUE,
                                                  width=12,
                                                  collapsible=TRUE,
                                                  {plotOutput('features_heatmap-heatmap') %>% withSpinner()}),
                                 features_selector=box(title='Select genes',
                                                       status='info',
                                                       solidHeader=TRUE,
                                                       width=4,
                                                       collapsible=TRUE,
                                                       features_heatmap.features_of_interest.tb,
                                                       features_heatmap.show_feature_names.checkbox,
                                                       features_heatmap.seurat_cluster_set.dropdown,
                                                       features_heatmap.use_average_value.checkbox))

  ### assemble tab content
  features_heatmap.content <-tabItem(tabName='features_heatmap-tab',
                                    h1('Visualise expression of multiple features in clusters'),
                                    fluidRow(project_name_text_box.ui(id='features_heatmap', width=12)),
                                    fluidRow(features_heatmap.boxes$heatmap_plot,
                                             features_heatmap.boxes$features_selector))

  ## submit/configure data tab
  submit_data.tab <- menuItem(text='Configure', tabName='submit_data_tab', icon=icon('cogs'), badgeLabel='!', badgeColor='yellow', selected=TRUE)

  ### define ui elements

  ### define layout boxes

  ### assemble tab content
  load_dataset.content <- tabItem(tabName='submit_data_tab',
                                  h1('Select a Seurat object'),
                                  available_seurats.ui(id='load_dataset'))

  ## menu tab hyperlinks
  email_me.tab <- menuItem(text='mail me', href='mailto:christopher.barrington@crick.ac.uk?subject=[seurat-vis] Hello there', icon=icon('comment-dots'), newtab=FALSE)
  bug_report.tab <- menuItem(text='report a bug', href='mailto:christopher.barrington@crick.ac.uk?subject=[seurat-vis] I found a bug!', icon=icon('bug'), newtab=FALSE)
  feature_request.tab <- menuItem(text='suggest a feature', href='mailto:christopher.barrington@crick.ac.uk?subject=[seurat-vis] It would be cool if...', icon=icon('lightbulb'), newtab=FALSE)
  github_link.tab <- menuItem(text=sprintf('GitHub (version: %s)', packageVersion('seuratvis')), href='github.com/ChristopherBarrington/seuratvis', icon=icon('code-branch'), newtab=FALSE)

  # header definition
  css <- 'table.dataTable tr.active td, table.dataTable td.active {background-color: #3C8DBC !important;}'
  dashboard_header <- dashboardHeader(title='seurat-vis')

  # dashboard body definition
  list(cell_filtering.content,
       gene_highlighting.content,
       features_heatmap.content,
       load_dataset.content) %>%
    do.call(what=tabItems) %>%
    dashboardBody(rclipboardSetup(), tags$head(tags$style(HTML(css)))) -> dashboard_body

  # sidebar definition
  list(cell_filtering.tab,
       gene_highlighting.tab,
       # discover_markers.tab,
       # features_heatmap.tab,
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
