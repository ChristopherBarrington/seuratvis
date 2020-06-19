
.ui_tab.feature_highlighting <- function() {
  bquote({
    # feature highlighting on map tab
    tab <- 'gene_highlighting'

    # define layout boxes
    boxes <- list(cluster_dim_plot=box(title='Clustered map',
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
                                    feature_picker.ui(id=tab),
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

    # assemble tab content
    content <- tabItem(tabName=NS(namespace=tab, id='tab'),
                       h1('Highlight cell features on the map'),
                       fluidRow(project_name_text_box.ui(id=tab, width=5),
                                picked_feature_and_description_text_box.ui(id=tab, width=7)),
                       fluidRow(number_of_reads_text_box.ui(id=tab, width=2),
                                number_of_clusters_text_box.ui(id=tab, width=2),
                                number_of_cells_text_box.ui(id=tab, width=2),
                                number_of_genes_in_assay_text_box.ui(id=tab, width=2),
                                number_of_reads_per_cell_text_box.ui(id=tab, width=2, f=median),
                                number_of_genes_per_cell_text_box.ui(id=tab, width=2, f=median)),
                       fluidRow(boxes$cluster_dim_plot,
                                boxes$gene_expression_map,
                                boxes$expression_per_cluster),
                       fluidRow(boxes$gene_selector,
                                boxes$plot_options))

    # assign variables in the parent environment
    menuItem(text='Highlight features', tabName=NS(namespace=tab, id='tab'), icon=icon('search')) %>% list() %>% append(x=menus) -> menus
    contents %<>% append(list(content))})
}
