
.ui_tab.gene_modules <- function() {
  bquote({
    # feature highlighting on map tab
    tab <- 'gene_modules'

    ## assemble tab content
    tabItem(tabName=NS(namespace=tab, id='tab'),
            h1('Gene modules'),
            fluidRow(project_name_text_box.ui(id=tab, width=8),
                     number_of_reads_text_box.ui(id=tab, width=2),
                     number_of_cells_text_box.ui(id=tab, width=2)),
            fluidRow(gene_module_score_plot.ui(id=tab))) -> content

    # assign variables in the parent environment
    menuItem(text='Gene modules', tabName=NS(namespace=tab, id='tab'), icon=icon('layer-group'), badgeLabel='!', badgeColor='red') %>% list() %>% append(x=menus) -> menus
    contents %<>% append(list(content))})
}
