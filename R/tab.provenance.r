#'
#' 
provenance.tab <- function() {
  bquote({

    tab <- 'provenance_tab'
    menuItem(text='Provenance', icon=icon('history'), tabName=tab) -> menu_item

    tabItem(tabName=tab,
            h1('View the functions used to create this Seurat object'),
            fluidRow(dataset_info_text_box.ui(id=NS(tab, 'project_name'), width=12)),
            ace_editor.ui(id=NS(tab, 'editor'))) -> content

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
      renderUI({provenace_picker.ui(id=tab, seurat=seurat)})  -> server_output$right_sidebar.data_opts
      renderUI({p('No options')}) -> server_output$right_sidebar.plotting_opts}})

  # call the modules for this tab
  provenace_picker <- callModule(module=provenace_picker.server, id='', seurat=seurat)

  callModule(module=dataset_info_text_box.project_name, id='project_name', seurat=seurat)
  callModule(module=ace_editor.server, id='editor', display_text=provenace_picker$script)
}
