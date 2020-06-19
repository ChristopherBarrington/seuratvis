
.ui_tab.findmarkers_table <- function() {
  bquote({
    # feature highlighting on map tab
    tab <- 'findmarkers_table'

    ## assemble tab content
    filter_FindMarkersResults.ui(id=tab) %>%
      append(list(tabName=NS(namespace=tab, id='tab'))) %>%
      do.call(what=tabItem) -> content

    # assign variables in the parent environment
    menuItem(text='Check FindMarkers results', tabName=NS(namespace=tab, id='tab'), icon=icon('binoculars'), badgeLabel='!', badgeColor='red') %>% list() %>% append(x=menus) -> menus
    contents %<>% append(list(content))})
}
