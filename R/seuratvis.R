#' launches the seuratvis app
#'
#' @return shiny application object
#'
#' @import shiny
#'
#' @examples
#' 
#' \dontrun{
#' seuratvis()}
#' 
#' @export
#' 
seuratvis <- function()
  shinyApp(ui=shinyAppUI, server=shinyAppServer)

shinyAppUI <- function(...) {
  # tab definitions
  ## initialise lists for tab menus and contents
  menus <- list()
  contents <- list()

  ## get the menu tabs and contents
  eval(preprocessing.tab())
  eval(highlight_features.tab())
  eval(cluster_classification.tab())
  eval(provenance.tab())
  eval(available_seurats.tab())
  # eval(contact_links.menu())

  # header dropdown definition
  contact_links.list() %>%
    dropdownMenu(.list=.,
                 type='notifications',
                 badgeStatus=NULL, headerText='', icon=icon('github')) -> contact_header_dropdown

  # header definition
  logo_lg <- htmltools::HTML("<p style='font-size:26px'>seurat<b>vis</b></p>")
  logo_sm <- htmltools::HTML("<p style='font-size:26px'>s<b>v</b></p>")
  dashboardHeaderPlus(
    title=tagList(span(class='logo-lg', logo_lg),
                  span(class='logo-mini', logo_sm)),
    enable_rightsidebar=TRUE,
    contact_header_dropdown) -> dashboard_header

  # dashboard body definition
  tags$head(tags$style(HTML(text='table.dataTable tr.active td, table.dataTable td.active {background-color: #3C8DBC !important;}'))) -> cssDT
  tags$style(type='text/css', '#provenance_tab-editor-ace_editor {height: calc(65vh) !important;} #cell_filtering_tab--ace_verbatim_text_output-ace_editor {height: 150px !important;') -> cssAce # apply this to class `ace_editor`
  tags$style(type='text/css', '.boy, .girl {font-size: x-large} .boy {color: #347DC1} .girl {color: #CC6594') -> cssSex

  append(contents,
         list()) %>%
    do.call(what=tabItems) %>%
    dashboardBody(rclipboardSetup(), cssDT, cssAce, shinyDashboardThemes(theme='grey_dark')) -> dashboard_body

  # sidebar definition
  menus %>% append(list(actionButton(inputId='clickme', label='', icon=icon('user-secret')))) %>%
    sidebarMenu(id='left_sidebar') %>%
    dashboardSidebar() -> left_dashboard_sidebar

  # right sidebar definition
  rightSidebar(title='Right Sidebar',
               background='dark',
               rightSidebarTabContent(id='data_opts', title='Options', icon='wrench', active=TRUE, uiOutput(outputId='right_sidebar.data_opts')),
               rightSidebarTabContent(id='plotting_opts', title='Plotting', icon='palette', active=FALSE, uiOutput(outputId='right_sidebar.plotting_opts')),
               rightSidebarTabContent(id='config_opts', title='Configure', icon='dna', active=FALSE, uiOutput(outputId='right_sidebar.config_opts'))) -> right_dashboard_sidebar

  # assemble the final UI
  list(header=dashboard_header,
       sidebar=left_dashboard_sidebar,
       rightsidebar=right_dashboard_sidebar,
       body=dashboard_body,
       title='seuratvis',
       skin='blue',
       useShinyjs()) %>%
    do.call(what=dashboardPagePlus)
}

shinyAppServer <- function(input, output, session) {
  # ###############################################################################################
  # call servers for modules
  ## Seurat object interaction modules
  available_seurats <- callModule(module=available_seurats_tab.server, id='configuration_tab', server_input=input, server_output=output, server_session=session)
  seurat <- callModule(module=process_seurat.server, id='process_seurat', server_input=input, server_output=output, server_session=session, available_seurats=available_seurats)

  ## load the filter_seurat module
  # dataset_filtering <- reactiveValues(foo='bar') #callModule(module=cell_filtering.server, id='seuratvis', seurat=seurat)

  ## load the servers for the analysis windows (menuItem or menuSubItem from the sidebar)
  callModule(module=cell_filtering_tab.server, id='cell_filtering_tab', server_input=input, server_output=output, server_session=session, seurat=seurat)
  callModule(module=dimensionality_tab.server, id='dimensionality_tab', server_input=input, server_output=output, server_session=session, seurat=seurat)
  callModule(module=highlight_feature_tab.server, id='highlight_feature_tab', server_input=input, server_output=output, server_session=session, seurat=seurat)
  callModule(module=highlight_feature_and_clusters_tab.server, id='highlight_feature_and_clusters_tab', server_input=input, server_output=output, server_session=session, seurat=seurat)
  callModule(module=highlight_multiple_features.server, id='highlight_multiple_features_tab', server_input=input, server_output=output, server_session=session, seurat=seurat)
  callModule(module=findmarkers_results_tab.server, id='findmarkers_results_tab', server_input=input, server_output=output, server_session=session, seurat=seurat)
  callModule(module=gene_module_score_in_clusters_tab.server, id='gene_module_score_in_clusters_tab', server_input=input, server_output=output, server_session=session, seurat=seurat)
  callModule(module=gene_module_scores_in_a_cluster_tab.server, id='gene_module_scores_in_a_cluster_tab', server_input=input, server_output=output, server_session=session, seurat=seurat)
  callModule(module=provenance_tab.server, id='provenance_tab', server_input=input, server_output=output, server_session=session, seurat=seurat)

observeEvent(input$clickme, {
  # seurat$object <- rnorm(1)
  # seurat$object <- get('human_CS12', envir=globalenv())
  # seurat$metadata <- seurat$object@meta.data
  # seurat$gene_module_scores <- select_at(seurat$metadata, vars(starts_with('GeneModule-')))
})

  # ###############################################################################################
  # reactions to tab selection
  ## react to opening tab with a filtered object loaded
  # observeEvent(input$sidebarmenu, {
  #   if(!is.null(seurat$object) & input$sidebarmenu=='cell_filtering-tab' && (!is.null(seurat$object@misc$cells_filtered) && seurat$object@misc$cells_filtered))
  #     sendSweetAlert(session=session, type='success', html=TRUE,
  #                    title='Notice', btn_labels='Great!',
  #                    text=tags$span('It looks like low-quality cells have already been removed from this Seurat object:', tags$h5(tags$code('@misc$cells_filtered == TRUE'))),
  #                    closeOnClickOutside=TRUE, showCloseButton=FALSE)})

  ## when a tab is selected from the left sidebar, activate the right sidebar and select the data_opts tab
  observeEvent(eventExpr=input$left_sidebar, handlerExpr={
    # open the sidebar
    removeClass(selector='body', class='control-sidebar-closed')
    addClass(selector='body', class='control-sidebar-open')

    # activate the data sidebar tab (or config sidebar tab if on the available seurats tab)
    # runjs("$('.nav-tabs a[href=\"#control-sidebar-plotting_opts-tab\"]').tab('show');")
    runjs("$('.nav-tabs a[href=\"#control-sidebar-data_opts-tab\"]').tab('show');")
    if(input$left_sidebar=='configuration_tab')
      runjs("$('.nav-tabs a[href=\"#control-sidebar-config_opts-tab\"]').tab('show');")

    if(input$left_sidebar!='configuration_tab' & !isTruthy(seurat$object)) {
      sendSweetAlert(session=session, title='SOS', btn_labels='OK', html=TRUE, closeOnClickOutside=TRUE, type='error',
                     text=tags$span('Select One Seurat', br(), br(), 'Use the config tab of the right sidebar to select a Seurat object to use in the app.'))
      updateNavbarPage(session=session, inputId='left_sidebar', selected='configuration_tab')
    }
  })

  # ###############################################################################################
  # any code to exectue when the session ends
  session$onSessionEnded(function() {
    message('### session ended')})
}

