
.ui_tab.configure_seurat <- function() {
  bquote({
    # submit/configure data tab
    tab <- 'configure_seurat'

    # assemble tab content
    content <- tabItem(tabName=NS(namespace=tab, id='tab'),
                       h1('Select a Seurat object'),
                       available_seurats.ui(id='load_dataset'))

    # assign variables in the parent environment
    menuItem(text='Configure', tabName=NS(namespace=tab, id='tab'), icon=icon('cogs'), selected=TRUE) %>% list() %>% append(x=menus) -> menus
    contents %<>% append(list(content))})
}
