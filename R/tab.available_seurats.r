#'
#' 
available_seurats.tab <- function() {
  bquote({
    menuItem(text='Configuration', icon=icon('dna'), startExpanded=TRUE,
             menuSubItem(text='Object information', tabName='configuration_tab', icon=menuSubItem_icon(), selected=TRUE),
             menuSubItem(text='Add feature modules', tabName='add_feature_module_tab', icon=menuSubItem_icon())) %>%
      modify_stop_propagation() -> menu_item

    list(tabItem(tabName='configuration_tab',
                 h1('Available Seurat objects'),
                 h5('These objects have been found in your workspace and can be loaded'),
                 fluidRow(boxPlus(title='', closable=FALSE, width=12, status='primary',
                                  seurats_in_workspace.table(id=NS('configuration_tab', 'seurats_table'))))),
         
         tabItem(tabName='add_feature_module_tab',
                 h1('Add feature modules'),
                 h5('Add a list of features and a feature module name below to add a feature module score to the Seurat object'),
                 add_feature_module_score.ui(id=NS('add_feature_module_tab', 'add_feature_module')))) -> content

    menus %<>% append(list(menu_item))
    contents %<>% append(content)})
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
  return(callModule(module=seurats_in_workspace.server, id='seurats_table'))
}

#'
#' 
add_feature_module_tab.server <- function(input, output, session, server_input, server_output, server_session, seurat) {
  # build the sidebar ui
  observeEvent(eventExpr=server_input$left_sidebar, handlerExpr={
    tab <- 'add_feature_module_tab'
    if(server_input$left_sidebar==tab) {
      tab %<>% str_c('-')
      renderUI({p('No options')})  -> server_output$right_sidebar.data_opts
      renderUI({p('No options')}) -> server_output$right_sidebar.plotting_opts}})

  # call the modules for this tab
  callModule(module=add_feature_module_score.server, id='add_feature_module', seurat=seurat)
}







