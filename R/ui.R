
shinyAppUI <- function(...) {
  # tab definitions
  ## initialise lists for tab menus and contents
  menus <- list()
  contents <- list()

  ## get the menu tabs and contents
  #! TODO: can this automatically get all `.ui_tab.` functions and run them?
  eval(.ui_tab.cell_filtering())
  eval(.ui_tab.feature_highlighting())
  eval(.ui_tab.feature_highlighting_cluster_selection())
  eval(.ui_tab.configure_seurat())
  eval(.ui_tab.sidebar_links())

  # header definition
  logo <- htmltools::HTML("<p style='font-size:26px'>seurat<b>vis</b></p>")
  dashboard_header <- dashboardHeader(title=logo)

  # dashboard body definition
  css <- 'table.dataTable tr.active td, table.dataTable td.active {background-color: #3C8DBC !important;}'
  append(contents,
         list()) %>%
    do.call(what=tabItems) %>%
    dashboardBody(rclipboardSetup(), tags$head(tags$style(HTML(css)))) -> dashboard_body

  # sidebar definition
  append(menus, 
         list()) %>%
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
