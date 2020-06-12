
.ui_tab.dimension_parameter_choice <- function() {
  bquote({
    # feature highlighting on map tab
    tab <- 'dimension_parameter_choice'

    # define layout boxes
    boxes <- list(elbow=box(title='Elbow plot',
                            footer='Variance described by successive principal components',
                            status='primary', solidHeader=TRUE, width=4, collapsible=TRUE,
                            elbow_plot.ui(id=tab)),
                  jackstraw=box(title='JackStraw plot',
                                footer='Compare random and observed PCA for features',
                                status='primary', solidHeader=TRUE, width=4, collapsible=TRUE,
                                jackstraw_plot.ui(id=tab)),
                  jackstraw_pvalues=box(title='JackStraw p-values',
                                        footer='Permutation test of component against random background',
                                        status='primary', solidHeader=TRUE, width=4, collapsible=TRUE,
                                        jackstraw_pvalue_plot.ui(id=tab)),
                  heatmap=box(title='Top contributing features',
                              footer='Contribution of features to a principle component',
                              status='primary', solidHeader=TRUE, width=4, collapsible=TRUE,
                              pca_heatmap.ui(id=tab)),
                  options=box(title='Plot options',
                              status='primary', solidHeader=TRUE, width=3, collapsible=TRUE,
                              reduction_method_picker.ui(id=tab, regex='^pca'),
                              principle_component_picker.ui(id=tab)))

    # assemble tab content
    content <- tabItem(tabName=NS(namespace=tab, id='tab'),
                       h1('Dataset dimensionality'),
                       project_name_text_box.ui(id=tab, width=12),
                       boxes$options,
                       column(width=8, fluidRow(boxes$elbow, boxes$heatmap),
                                       fluidRow(boxes$jackstraw, boxes$jackstraw_pvalues)))

    content <- tabItem(tabName=NS(namespace=tab, id='tab'),
                       h1('Dataset dimensionality'),
                       fluidRow(project_name_text_box.ui(id=tab, width=12)),
                       fluidRow(boxes$elbow, boxes$heatmap, boxes$options),
                       fluidRow(boxes$jackstraw, boxes$jackstraw_pvalues))

    # assign variables in the parent environment
    menuItem(text='Dataset dimensionality', tabName=NS(namespace=tab, id='tab'), icon=icon('bar-chart')) %>% list() %>% append(x=menus) -> menus
    contents %<>% append(list(content))})
}
