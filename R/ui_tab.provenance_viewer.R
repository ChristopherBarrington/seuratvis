
.ui_tab.provenance_viewer <- function() {
  bquote({
    # provenance viewer tab
    tab <- 'provenance_viewer'

    # define ui-specific css
    cssAce <- tags$style(type='text/css', '#r_script {height: calc(100vh-80px) !important;}')


    # assemble tab content
    content <- tabItem(tabName=NS(namespace=tab, id='tab'),
                       h1('View !the functions used to create this Seurat object'),
                       fluidRow(project_name_text_box.ui(id=tab, width=6),
                                number_of_reads_text_box.ui(id=tab, width=3),
                                number_of_cells_text_box.ui(id=tab, width=3)),
                       provenance_step_viewer.ui(id=tab))

    # assign variables in the parent environment
    menuItem(text='Provenance', tabName=NS(namespace=tab, id='tab'), icon=icon('history')) %>% list() %>% append(x=menus) -> menus
    contents %<>% append(list(content))})
}
