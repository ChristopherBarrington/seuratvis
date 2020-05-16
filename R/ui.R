
# library(shiny)
# library(shinydashboard)
# library(dqshiny)
# library(shinycssloaders)
# library(ggvis)
# library(tidyverse)
# library(magrittr)

# load('int.RData') ######
# # seurat <- human_CS17_thoracic
# seurat <- human_CS17_brachial
# starter_gene <- sample(x=rownames(seurat), size=1)
# starter_gene <- 'SOX2'
# ################

# tab definitions

## cell filtering tab
cell_filtering.tab <- menuItem(text='Cell filtering', tabName='cell_filtering-tab', icon=icon('filter'))

### define ui elements
min_features_per_cell.textinput <- textInput(inputId='min_features_per_cell.textinput', label='Minimum features per cell', placeholder='min')
max_features_per_cell.textinput <- textInput(inputId='max_features_per_cell.textinput', label='Maximum features per cell', placeholder='max')

min_expression_per_cell.textinput <- textInput(inputId='min_expression_per_cell.textinput', label='Minimum total UMIs per cell', placeholder='min')
max_expression_per_cell.textinput <- textInput(inputId='max_expression_per_cell.textinput', label='Maximum total UMIs per cell', placeholder='max')

percent_mitochondria.textinput <- textInput(inputId='percent_mitochondria.textinput', label='Maximum mitochondrial (%)', placeholder='max')

subset_conditions.textoutput <- verbatimTextOutput(outputId='cell_filtering-subset_conditions')
subset_conditions.plain.copybutton <- uiOutput(outputId='cell_filtering-subset_conditions.plain', inline=TRUE)
subset_conditions.tsv.copybutton <- uiOutput(outputId='cell_filtering-subset_conditions.tsv', inline=TRUE)
subset_conditions.r.copybutton <- uiOutput(outputId='cell_filtering-subset_conditions.r', inline=TRUE)

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
                                                    column(width=2, min_expression_per_cell.textinput, max_expression_per_cell.textinput),
                                                    column(width=2, min_features_per_cell.textinput, max_features_per_cell.textinput),
                                                    column(width=2, percent_mitochondria.textinput),
                                                    column(width=4, shiny::tags$label('Conditional expression to select cells'), br(), subset_conditions.textoutput,
                                                                    shiny::tags$label('Copy subset conditions to clipboard'), br(), subset_conditions.plain.copybutton, subset_conditions.tsv.copybutton, subset_conditions.r.copybutton))) %>%
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
gene_highlighting.tab <- menuItem(text='Highlight genes', tabName='gene_highlighting-tab', icon=icon('search'))

### define ui elements
gene_name.dropdown <- autocomplete_input(id='gene_of_interest.dd', label='Gene name',
                                         # options=sort(rownames(seurat)),
                                         options=NULL,
                                         placeholder='Gene name', value='orig.ident')

seurat_cluster_set.dropdown <- selectInput(inputId='seurat_cluster_set.dd', label='Seurat cluster definitions',
                                           choices='seurat_clusters', selected='seurat_clusters',
                                           multiple=FALSE)
gene_highlighting.reduction_selection.dd <- selectInput(inputId='reduction_selection.dd', label='Reduction method',
                                                        choices=NULL, selected=NULL,
                                                        multiple=FALSE)
gene_highlighting.assay_selection.dd <- selectInput(inputId='assay_selection.dd', label='Assay',
                                                    choices=NULL, selected=NULL,
                                                    multiple=FALSE)
expression_range.slider <- sliderInput(inputId='expression_range.slider', label='Expression limits',
                                       # min=0, max=round(max(FetchData(seurat, starter_gene))+0.05), step=0.1, value=c(-Inf,Inf))
                                       min=0, max=1, step=0.1, value=c(-Inf,Inf))

opacity.slider <- sliderInput(inputId='opacity.slider', label='Opacity', min=0.1, max=1, step=0.1, value=1)

gene_highlighting.point_size <- sliderInput(inputId='gene_highlighting.point_size.slider', label='Point size',
                                            min=0.1, max=3, step=0.1, value=0.6)
                                            # min=1, max=15, step=1, value=3)

gene_highlighting.label_clusters <- checkboxInput(inputId='gene_highlighting.label_clusters.checkbox', label='Label clusters', value=TRUE)
gene_highlighting.label_clusters <- materialSwitch(inputId='gene_highlighting.label_clusters.checkbox', label='Cluster labels', value=TRUE, right=TRUE, status='success')

#### colour selector palette box
# colourInput.defaults <- list(showColour='both', palette='limited', allowedCols=default_colour_palette(), allowTransparent=FALSE, returnName=TRUE)
# append(colourInput.defaults, list(inputId='expression_min.colour', label='Low', value='linen')) %>% do.call(what=colourInput) -> expression_min.colour.selector
# append(colourInput.defaults, list(inputId='expression_max.colour', label='High', value='darkviolet')) %>% do.call(what=colourInput) -> expression_max.colour.selector
# materialSwitch(inputId='expression_palette_full', label='Show full palette?', value=FALSE, right=TRUE, status='success') -> expression_paletette_type_selector

### define layout boxes
gene_highlighting.boxes <- list(cluster_dim_plot=box(title='Clustered map',
                                                     footer='Seurat map showing cell clusters',
                                                     status='success',
                                                     solidHeader=TRUE,
                                                     width=4,
                                                     collapsible=TRUE,
                                                     {plotOutput('genes_highlighting-clustered_map') %>% withSpinner()}),
                                gene_expression_map=box(title='Marker gene on map',
                                                        footer='Cells coloured by expression',
                                                        status='success',
                                                        solidHeader=TRUE,
                                                        width=4,
                                                        collapsible=TRUE,
                                                        {plotOutput('genes_highlighting-gene_expression_map') %>% withSpinner()}),
                                expression_per_cluster=box(title='Marker gene in clusters',
                                                           footer='Expression of gene in clusters',
                                                           status='success',
                                                           solidHeader=TRUE,
                                                           width=4,
                                                           collapsible=TRUE,
                                                           {plotOutput('genes_highlighting-expression_per_cluster') %>% withSpinner()}),
                                                           # {plotlyOutput('genes_highlighting-expression_per_cluster') %>% withSpinner()}),
                                gene_selector=box(title='Select gene',
                                                  status='primary',
                                                  solidHeader=TRUE,
                                                  width=4,
                                                  collapsible=TRUE,
                                                  gene_name.dropdown,
                                                  expression_range.slider,
                                                  seurat_cluster_set.dropdown,
                                                  gene_highlighting.label_clusters,
                                                  gene_highlighting.reduction_selection.dd,
                                                  gene_highlighting.assay_selection.dd),
                                plot_options=box(title='Plot options',
                                                 status='primary',
                                                 solidHeader=TRUE,
                                                 width=4,
                                                 collapsible=TRUE,
                                                 shiny::tags$label('Feature value colours'), br(),
                                                 colour_palette.ui(id='gene_highlighting', include_full=TRUE,
                                                                   selectors=list(list(inputId='low', label='Low', value='linen'),
                                                                                  list(inputId='high', label='High', value='darkviolet'))),
                                                 gene_highlighting.point_size,
                                                 opacity.slider))

### assemble tab content
gene_highlighting.content <- tabItem(tabName='gene_highlighting-tab',
                                     h1('Highlight expression of interesting genes on the map'),
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
prettyRadioButtons(inputId='seurat_select.input', label='Select a Seurat object', 
                   choiceNames='', choiceValues='',
                   icon=icon('check'), shape='curve', outline=TRUE, bigger=TRUE, status='primary', animation='smooth') -> seurat_select.checkbox

### define layout boxes

### assemble tab content
load_dataset.content <- tabItem(tabName='submit_data_tab',
                                h1('Select a Seurat object'),
                                seurat_select.checkbox)

## menu tab hyperlinks
email_me.tab <- menuItem(text='mail me', href='mailto:christopher.barrington@crick.ac.uk?subject=[seurat-vis] Hello there', icon=icon('comment-dots'))
bug_report.tab <- menuItem(text='report a bug', href='mailto:christopher.barrington@crick.ac.uk?subject=[seurat-vis] I found a bug!', icon=icon('bug'))
feature_request.tab <- menuItem(text='suggest a feature', href='mailto:christopher.barrington@crick.ac.uk?subject=[seurat-vis] It would be cool if...', icon=icon('lightbulb'))
github_link.tab <- menuItem(text=sprintf('GitHub (version: %s)', packageVersion('seuratvis')), href='github.com/ChristopherBarrington/seuratvis', icon=icon('code-branch'))

# header definition
dashboard_header <- dashboardHeader(title='seurat-vis')

# dashboard body definition
list(cell_filtering.content,
     gene_highlighting.content,
     features_heatmap.content,
     load_dataset.content) %>%
  do.call(what=tabItems) %>%
  dashboardBody(rclipboardSetup()) -> dashboard_body

# sidebar definition
list(cell_filtering.tab,
     gene_highlighting.tab,
     discover_markers.tab,
     features_heatmap.tab,
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
  do.call(what=dashboardPage) -> shinyAppUI
