#'
#' 
configuration.tab <- function() {
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
configuration_tab.server <- function(input, output, session, server_input, server_output, server_session) {
  # build the sidebar ui
  observeEvent(eventExpr=server_input$left_sidebar, handlerExpr={
    tab <- 'configuration_tab'
    
    if(server_input$left_sidebar==tab) {
      tab %<>% str_c('-')
      addClass(selector='body', class='control-sidebar-closed')
      renderUI({p('No options')})  -> server_output$right_sidebar.data_opts
      renderUI({p('No options')}) -> server_output$right_sidebar.plotting_opts
      renderUI({p('No options')}) -> server_output$right_sidebar.config_opts}})

  # call the modules for this tab
  # provenace_picker <- callModule(module=provenace_picker.server, id='', seurat=seurat)

  # callModule(module=project_name_text_box.server, id='project_name', seurat=seurat)
  # callModule(module=provenance_step_viewer.server, id='editor', picked_provenance=provenace_picker)
}

seurats_in_workspace.table <- function(...) return(h2('A TABLE'))
