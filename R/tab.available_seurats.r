#'
#' 
available_seurats.tab <- function() {
  bquote({

    tab <- 'configuration_tab'
    menuItem(text='Configure', tabName=tab, icon=icon('map'), selected=TRUE) -> menu_item

    tabItem(tabName='configuration_tab',
            h1('Available Seurat objects'),
            h5('These objects have been found in your workspace and can be loaded'),
            fluidRow(boxPlus(title='', closable=FALSE, width=12, status='primary',
                             seurats_in_workspace.table(id=NS(tab, 'seurats_table'))))) -> content

    menus %<>% append(list(menu_item))
    contents %<>% append(list(content))})
}

#'
#' 
available_seurats_tab.server <- function(input, output, session, server_input, server_output, server_session) {
  # build the sidebar ui
  observeEvent(eventExpr=server_input$left_sidebar, handlerExpr={
    tab <- 'configuration_tab'
    if(server_input$left_sidebar==tab) {
      tab %<>% str_c('-')
      renderUI({p('No options')})  -> server_output$right_sidebar.data_opts
      renderUI({p('No options')}) -> server_output$right_sidebar.plotting_opts}})

  # call the modules for this tab
  return(callModule(seurats_in_workspace.server, id='seurats_table'))
}

