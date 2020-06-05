
.ui_tab.cell_filtering <- function(...) {
  bquote({
    # cell filtering tab
    tab <- 'cell_filtering'

    # define layout boxes
    plot_boxes.defaults <- list(solidHeader=TRUE, width=12, collapsible=TRUE, collapsed=FALSE)
    boxes <- list(total_expression=append(plot_boxes.defaults,
                                          list(title='UMI per cell', footer='Total UMI attributed to a cell', status='primary',
                                               column(width=2, offset=1, boxplot_plot.ui(id=tab, feature='nCount_RNA')),
                                               column(width=3, offset=0, knee_plot.ui(id=tab, feature='nCount_RNA')),
                                               column(width=4, offset=1, density_plot.ui(id=tab, feature='nCount_RNA')))) %>%
                                   do.call(what=box),

                  unique_features=append(plot_boxes.defaults,
                                         list(title='Features per cell', footer='Number of distinct features detected in a cell', status='primary',
                                              column(width=2, offset=1, boxplot_plot.ui(id=tab, feature='nFeature_RNA')),
                                              column(width=3, offset=0, knee_plot.ui(id=tab, feature='nFeature_RNA')),
                                              column(width=4, offset=1, density_plot.ui(id=tab, feature='nFeature_RNA')))) %>%
                                  do.call(what=box),

                  percent_mitochondria=append(plot_boxes.defaults,
                                              list(title='Mitochondrial UMI', footer='Proportion of mitochondrial genes detected in a cell', status='primary',
                                                   column(width=2, offset=1, boxplot_plot.ui(id=tab, feature='percent_mt')),
                                                   column(width=3, offset=0, knee_plot.ui(id=tab, feature='percent_mt')),
                                                   column(width=4, offset=1, density_plot.ui(id=tab, feature='percent_mt')))) %>%
                                       do.call(what=box),

                  thresholds=append(plot_boxes.defaults,
                                    list(title='Thresholds', status='success',
                                         column(width=2, total_umi_per_cell_filter.ui(id=tab)),
                                         column(width=2, features_per_cell_filter.ui(id=tab)),
                                         column(width=2, percent_mt_per_cell_filter.ui(id=tab)),
                                         column(width=4, show_filtering_parameters.ui(id=tab)))) %>%
                             do.call(what=box))

    # assemble tab content
    content <- tabItem(tabName=NS(namespace=tab, id='tab'),
                       h1('Identify cells that can be removed with quality filters'),
                       fluidRow(project_name_text_box.ui(id=tab, width=4),
                                number_of_reads_text_box.ui(id=tab, width=2),
                                number_of_cells_text_box.ui(id=tab, width=2),
                                number_of_reads_per_cell_text_box.ui(id=tab, width=2, f=median),
                                number_of_genes_per_cell_text_box.ui(id=tab, width=2, f=median)),
                       fluidRow(boxes$total_expression),
                       fluidRow(boxes$unique_features),
                       fluidRow(boxes$percent_mitochondria),
                       fluidRow(boxes$thresholds))

    # assign variables in the parent environment
    menuItem(text='Cell filtering', tabName=NS(namespace=tab, id='tab'), icon=icon('filter')) %>% list() %>% append(x=menus) -> menus
    contents %<>% append(list(content))})
}
