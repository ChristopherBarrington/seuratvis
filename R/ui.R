
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
starter_gene <- sample(x=rownames(seurat), size=1)

# ################


# tab definitions

## cell filtering tab
cell_filtering.tab <- menuItem(text='Cell filtering', tabName='cell_filtering-tab', icon=icon('filter'))

### define ui elements
percent_mitochondria.slider <- sliderInput(inputId='percent_mitochondria.slider', label='Mitochondria fraction',
                                           min=0,
                                           max=round(max(FetchData(seurat, 'percent_mt'))+(0.1/2)),
                                           value=round(max(FetchData(seurat, 'percent_mt'))+(0.1/2)),
                                           step=0.1)

step_size <- 100
min_genes_per_cell.slider <- sliderInput(inputId='min_genes_per_cell.slider', label='Minimum number of genes per cell',
                                         min={FetchData(seurat, 'nFeature_RNA') %>% min() %>% divide_by(step_size) %>% floor() %>% multiply_by(step_size)},
                                         max={FetchData(seurat, 'nFeature_RNA') %>% deframe() %>% quantile(p=0.25) %>% divide_by(step_size) %>% ceiling() %>% multiply_by(step_size)},
                                         value={FetchData(seurat, 'nFeature_RNA') %>% min() %>% divide_by(step_size) %>% floor() %>% multiply_by(step_size)},
                                         step=step_size)

max_genes_per_cell.slider <- sliderInput(inputId='max_genes_per_cell.slider', label='Maximum number of genes per cell',
                                         min={FetchData(seurat, 'nFeature_RNA') %>% deframe() %>% quantile(p=0.85) %>% divide_by(step_size) %>% ceiling() %>% multiply_by(step_size)},
                                         max={FetchData(seurat, 'nFeature_RNA') %>% max() %>% divide_by(step_size) %>% ceiling() %>% multiply_by(step_size)},
                                         value={FetchData(seurat, 'nFeature_RNA') %>% max() %>% divide_by(step_size) %>% ceiling() %>% multiply_by(step_size)},
                                         step=step_size)

min_expression_per_cell.slider <- sliderInput(inputId='min_expression_per_cell.slider', label='Minimum total expression per cell',
                                              min={FetchData(seurat, 'nCount_RNA') %>% min() %>% divide_by(step_size) %>% floor() %>% multiply_by(step_size)},
                                              max={FetchData(seurat, 'nCount_RNA') %>% deframe() %>% quantile(p=0.25) %>% divide_by(step_size) %>% ceiling() %>% multiply_by(step_size)},
                                              value={FetchData(seurat, 'nCount_RNA') %>% min() %>% divide_by(step_size) %>% floor() %>% multiply_by(step_size)},
                                              step=step_size)

max_expression_per_cell.slider <- sliderInput(inputId='max_expression_per_cell.slider', label='Maximum total expression per cell',
                                              min={FetchData(seurat, 'nCount_RNA') %>% deframe() %>% quantile(p=0.85) %>% divide_by(step_size) %>% ceiling() %>% multiply_by(step_size)},
                                              max={FetchData(seurat, 'nCount_RNA') %>% max() %>% divide_by(step_size) %>% ceiling() %>% multiply_by(step_size)},
                                              value={FetchData(seurat, 'nCount_RNA') %>% max() %>% divide_by(step_size) %>% ceiling() %>% multiply_by(step_size)},
                                              step=step_size)

### define layout boxes
cell_filtering.plot_boxes.defaults <- list(status='success', solidHeader=TRUE, width=4, collapsible=TRUE)
cell_filtering.plot_boxes.total_expression.defaults <- append(cell_filtering.plot_boxes.defaults, list(title='Expression per cell', footer='Total number of reads attributed to a cell'))
cell_filtering.plot_boxes.unique_genes.defaults <- append(cell_filtering.plot_boxes.defaults, list(title='Genes per cell', footer='Number of distinct genes detected in a cell'))
cell_filtering.plot_boxes.percent_mitochondria.defaults <- append(cell_filtering.plot_boxes.defaults, list(title='Mitochondrial expression', footer='Proportion of mitochondrial genes detected in a cell'))
cell_filtering.boxes <- list(total_expression_knee={append(cell_filtering.plot_boxes.total_expression.defaults, list({plotOutput(outputId='cell_filtering-total_expression_knee') %>% withSpinner()})) %>% do.call(what=box)},
                             total_expression_density={append(cell_filtering.plot_boxes.total_expression.defaults, list({plotOutput(outputId='cell_filtering-total_expression_density') %>% withSpinner()})) %>% do.call(what=box)},
                             total_expression_boxplot={append(cell_filtering.plot_boxes.total_expression.defaults, list({plotOutput(outputId='cell_filtering-total_expression_boxplot') %>% withSpinner()})) %>% do.call(what=box)},

                             unique_features_knee={append(cell_filtering.plot_boxes.unique_genes.defaults, list({plotOutput(outputId='cell_filtering-unique_genes_knee') %>% withSpinner()})) %>% do.call(what=box)},
                             unique_features_density={append(cell_filtering.plot_boxes.unique_genes.defaults, list({plotOutput(outputId='cell_filtering-unique_genes_density') %>% withSpinner()})) %>% do.call(what=box)},
                             unique_features_boxplot={append(cell_filtering.plot_boxes.unique_genes.defaults, list({plotOutput(outputId='cell_filtering-unique_genes_boxplot') %>% withSpinner()})) %>% do.call(what=box)},

                             percent_mitochondria_knee={append(cell_filtering.plot_boxes.percent_mitochondria.defaults, list({plotOutput(outputId='cell_filtering-percent_mitochondria_knee') %>% withSpinner()})) %>% do.call(what=box)},
                             percent_mitochondria_density={append(cell_filtering.plot_boxes.percent_mitochondria.defaults, list({plotOutput(outputId='cell_filtering-percent_mitochondria_density') %>% withSpinner()})) %>% do.call(what=box)},
                             percent_mitochondria_boxplot={append(cell_filtering.plot_boxes.percent_mitochondria.defaults, list({plotOutput(outputId='cell_filtering-percent_mitochondria_boxplot') %>% withSpinner()})) %>% do.call(what=box)},

                             thresholds=box(title='Thresholds',
                                            status='info',
                                            solidHeader=TRUE,
                                            width=12,
                                            collapsible=TRUE,
                                            column(width=4, min_expression_per_cell.slider, max_expression_per_cell.slider),
                                            column(width=4, min_genes_per_cell.slider, max_genes_per_cell.slider),
                                            column(width=4, percent_mitochondria.slider)),
                             project_name=valueBoxOutput(outputId='cell_filtering.project_name_box', width=4),
                             n_reads=valueBoxOutput(outputId='cell_filtering.n_reads_box', width=2),
                             n_cells=valueBoxOutput(outputId='cell_filtering.n_cells_box', width=2),
                             n_genes_per_cell=valueBoxOutput(outputId='cell_filtering.n_genes_per_cell_box', width=2),
                             n_reads_per_cell=valueBoxOutput(outputId='cell_filtering.n_reads_per_cell_box', width=2))

### assemble tab content
# cell_filtering.content <- tabItem(tabName='cell_filtering-tab',
#                                   h1('Identify cells that can be removed with quality filters'),
#                                   fluidRow(cell_filtering.boxes$project_name,
#                                            cell_filtering.boxes$n_reads,
#                                            cell_filtering.boxes$n_cells,
#                                            cell_filtering.boxes$n_reads_per_cell,
#                                            cell_filtering.boxes$n_genes_per_cell))
cell_filtering.content <- tabItem(tabName='cell_filtering-tab',
                                  h1('Identify cells that can be removed with quality filters'),
                                  fluidRow(cell_filtering.boxes$project_name,
                                           cell_filtering.boxes$n_reads,
                                           cell_filtering.boxes$n_cells,
                                           cell_filtering.boxes$n_reads_per_cell,
                                           cell_filtering.boxes$n_genes_per_cell),
                                  fluidRow(cell_filtering.boxes$total_expression_boxplot,
                                           cell_filtering.boxes$total_expression_knee,
                                           cell_filtering.boxes$total_expression_density),
                                  fluidRow(cell_filtering.boxes$unique_features_boxplot,
                                           cell_filtering.boxes$unique_features_knee,
                                           cell_filtering.boxes$unique_features_density),
                                  fluidRow(cell_filtering.boxes$percent_mitochondria_boxplot,
                                           cell_filtering.boxes$percent_mitochondria_knee,
                                           cell_filtering.boxes$percent_mitochondria_density),
                                  fluidRow(cell_filtering.boxes$thresholds))


## gene highlighting on map tab
gene_highlighting.tab <- menuItem(text='Highlight genes', tabName='gene_highlighting-tab', icon=icon('search'))

### define ui elements
gene_name.dropdown <- autocomplete_input(id = 'gene_of_interest.dd', label = 'Gene name',
                                         options = sort(rownames(seurat)),
                                         placeholder = 'Gene name', value = starter_gene)

seurat_cluster_set.dropdown <- selectInput(inputId = 'seurat_cluster_set.dd', label = 'Seurat cluster definitions',
                                           choices = 'seurat_clusters', selected = 'seurat_clusters',
                                           multiple = FALSE)
#### colour selector palette box
colourInput.defaults <- list(showColour='both', palette='limited', allowedCols=colour_palette, allowTransparent=FALSE, returnName=TRUE)
append(colourInput.defaults, list(inputId='expression_min.colour', label='Low', value='linen')) %>% do.call(what=colourInput) -> expression_min.colour.selector
append(colourInput.defaults, list(inputId='expression_max.colour', label='High', value='darkviolet')) %>% do.call(what=colourInput) -> expression_max.colour.selector
materialSwitch(inputId='expression_palette_full', label='Show full palette?', value=FALSE, right=TRUE, status='success') -> expression_paletette_type_selector

expression_range.slider <- sliderInput(inputId='expression_range.slider', label='Expression limits',
                                       min=0, max=round(max(FetchData(seurat, starter_gene))+0.05), step=0.1, value=c(-Inf,Inf))

opacity.slider <- sliderInput(inputId='opacity.slider', label='Opacity', min=0.1, max=1, step=0.1, value=1)

gene_highlighting.point_size <- sliderInput(inputId='gene_highlighting.point_size.slider', label='Point size',
                                            min=0.1, max=3, step=0.1, value=0.6)
                                            # min=1, max=15, step=1, value=3)

gene_highlighting.label_clusters <- checkboxInput(inputId='gene_highlighting.label_clusters.checkbox', label='Label clusters', value=TRUE)

### define layout boxes
gene_highlighting.boxes <- list(cluster_dim_plot=box(title='Clustered map',
                                                     footer='Seurat map showing cell clusters',
                                                     status='success',
                                                     solidHeader=TRUE,
                                                     width=4,
                                                     collapsible=TRUE,
                                                     # ggvisOutput('genes_highlighting-clustered_map'),
                                                     {plotOutput('genes_highlighting-clustered_map') %>% withSpinner()}),
                                gene_expression_map=box(title='Marker gene on map',
                                                        footer='Cells coloured by expression',
                                                        status='success',
                                                        solidHeader=TRUE,
                                                        width=4,
                                                        collapsible=TRUE,
                                                        # ggvisOutput('genes_highlighting-gene_expression_map'),
                                                        {plotOutput('genes_highlighting-gene_expression_map') %>% withSpinner()}),
                                expression_per_cluster=box(title='Marker gene in clusters',
                                                           footer='Expression of gene in clusters',
                                                           status='success',
                                                           solidHeader=TRUE,
                                                           width=4,
                                                           collapsible=TRUE,
                                                           # ggvisOutput('genes_highlighting-expression_per_cluster')),
                                                           {plotOutput('genes_highlighting-expression_per_cluster') %>% withSpinner()}),
                                gene_selector=box(title='Select gene',
                                                  status='info',
                                                  solidHeader=TRUE,
                                                  width=4,
                                                  collapsible=TRUE,
                                                  gene_name.dropdown,
                                                  expression_range.slider,
                                                  seurat_cluster_set.dropdown,
                                                  gene_highlighting.label_clusters),
                                plot_options=box(title='Plot options',
                                                 status='info',
                                                 solidHeader=TRUE,
                                                 width=4,
                                                 collapsible=TRUE,
                                                 shiny::tags$label('Feature value colours'), br(),
                                                 column(width=6, expression_min.colour.selector),
                                                 column(width=6, expression_max.colour.selector),
                                                 column(width=12, expression_paletette_type_selector),
                                                 gene_highlighting.point_size,
                                                 opacity.slider),
                                n_cells=valueBoxOutput(outputId='genes_highlighting.n_cells_box', width=2),
                                n_genes=valueBoxOutput(outputId='genes_highlighting.n_genes_box', width=2),
                                n_reads=valueBoxOutput(outputId='genes_highlighting.n_reads_box', width=2),
                                n_genes_per_cell=valueBoxOutput(outputId='genes_highlighting.n_genes_per_cell_box', width=2),
                                n_reads_per_cell=valueBoxOutput(outputId='genes_highlighting.n_reads_per_cell_box', width=2),
                                n_clusters=valueBoxOutput(outputId='genes_highlighting.n_clusters_box', width=2),
                                project_name=valueBoxOutput(outputId='genes_highlighting.project_name_box', width=5),
                                selected_gene=valueBoxOutput(outputId='genes_highlighting.selected_gene_box', width=7))

### assemble tab content
gene_highlighting.content <- tabItem(tabName='gene_highlighting-tab',
                                     h1('Highlight expression of interesting genes on the map'),
                                     fluidRow(gene_highlighting.boxes$project_name,
                                              gene_highlighting.boxes$selected_gene),
                                     fluidRow(gene_highlighting.boxes$n_reads,
                                              gene_highlighting.boxes$n_clusters,
                                              gene_highlighting.boxes$n_cells,
                                              gene_highlighting.boxes$n_genes,
                                              gene_highlighting.boxes$n_reads_per_cell,
                                              gene_highlighting.boxes$n_genes_per_cell),
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
features_heatmap.tab <- menuItem(text='Features heatmap', tabName='features_heatmap-tab', icon=icon('crosshairs'))

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
                                                     features_heatmap.use_average_value.checkbox),
                               project_name=valueBoxOutput(outputId='features_heatmap.project_name_box', width=12))

### assemble tab content
features_heatmap.content <-tabItem(tabName='features_heatmap-tab',
                                  h1('Visualise expression of multiple features in clusters'),
                                  fluidRow(features_heatmap.boxes$project_name),
                                  fluidRow(features_heatmap.boxes$heatmap_plot,
                                           features_heatmap.boxes$features_selector))

## submit/configure data tab
submit_data.tab <- menuItem(text='Submit data', tabName='submit_data_tab', icon=icon('cloud-upload'), badgeLabel='!', badgeColor='yellow', selected=TRUE)

### define ui elements

### define layout boxes

### assemble tab content
load_dataset.content <- tabItem(tabName='submit_data_tab',
                                h1('Load a Seurat object'))

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
  dashboardBody() -> dashboard_body

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
  sidebarMenu() %>%
  dashboardSidebar() -> dashboard_sidebar

# assemble the final UI
list(header=dashboard_header,
     sidebar=dashboard_sidebar,
     body=dashboard_body,
     title='seurat-vis',
     skin='blue') %>%
  do.call(what=dashboardPage) -> shinyAppUI

