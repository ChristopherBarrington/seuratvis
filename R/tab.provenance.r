#'
#' 
provenance.tab <- function() {
  bquote({

    tab <- 'provenance_tab'
    menuItem(text='Provenance', icon=icon('history'), tabName=tab) -> menu_item

    tabItem(tabName=tab,
            h1('View the functions used to create this Seurat object'),
            fluidRow(project_name_text_box.ui(id=NS(tab, 'project_name'), width=12)),
            provenance_step_viewer.ui(id=NS(tab, 'editor'))) -> content

    menus %<>% append(list(menu_item))
    contents %<>% append(list(content))})
}

#'
#' 
provenance_tab.server <- function(input, output, session, server_input, server_output, server_session, seurat) {
  # build the sidebar ui
  observeEvent(eventExpr=server_input$left_sidebar, handlerExpr={
    tab <- 'provenance_tab'
    
    if(server_input$left_sidebar==tab) {
      tab %<>% str_c('-')
      shinyjs::addClass(selector='body', class='control-sidebar-open')
      renderUI({provenace_picker.ui(id=tab, seurat=seurat)})  -> server_output$right_sidebar.data_opts
      renderUI({p('No options')}) -> server_output$right_sidebar.plotting_opts}})

  # call the modules for this tab
  provenace_picker <- callModule(module=provenace_picker.server, id='', seurat=seurat)

  callModule(module=project_name_text_box.server, id='project_name', seurat=seurat)
  callModule(module=provenance_step_viewer.server, id='editor', picked_provenance=provenace_picker)
}
