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
  eval(configuration.tab())
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
  tags$style(type='text/css', '.ace_editor {height: calc(65vh) !important;}') -> cssAce # apply this to class `ace_editor`
  tags$style(type='text/css', '.boy, .girl {font-size: x-large} .boy {color: #347DC1} .girl {color: #CC6594') -> cssSex

  append(contents,
         list()) %>%
    do.call(what=tabItems) %>%
    dashboardBody(cssDT, cssAce, shinyDashboardThemes(theme='grey_dark')) -> dashboard_body

  # sidebar definition
  menus %>%
    sidebarMenu(id='left_sidebar') %>%
    dashboardSidebar() -> left_dashboard_sidebar

  # right sidebar definition
  rightSidebar(title='Right Sidebar',
               background='dark',
               rightSidebarTabContent(id='data_opts', title='Options', icon='wrench', active=TRUE, uiOutput(outputId='right_sidebar.data_opts')),
               rightSidebarTabContent(id='plotting_opts', title='Plotting', icon='glasses', active=FALSE, uiOutput(outputId='right_sidebar.plotting_opts'))) -> right_dashboard_sidebar

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
  # scour the session for Seurat objects and populate the UI --------------------------------------
  #! TODO: this is not the right place for this. maybe a call to utils::globalVariables() in seuratvis() ... ?
  # available_seurat_objects <- find_seurat_objects() # has to be here to access the global environment ... ?

  # ###############################################################################################
  # call servers for modules
  ## Seurat object interaction modules
  #! TODO: move the find_seurats function to this module, and return the values as a rective. move to the seratvis id?
s <- get('human_CS12', envir=globalenv())
reactiveValues(
    object=s,
    formatted_project_name=Project(s) %>% reformat_project_name(),
    metadata=s@meta.data,
    features_in_assays=list(),
    reductions=Reductions(s),
    assays=Assays(s),
    gene_module_scores=select_at(s@meta.data, vars(starts_with('GeneModule-'))),
    gene_modules=s@misc$gene_modules,
    cluster_resolutions={resolutions <- c('seurat_clusters', str_subset(colnames(s@meta.data), '_snn_res.'))},
    all_idents={resolutions <- c('seurat_clusters', str_subset(colnames(s@meta.data), '_snn_res.')) ; select_at(s@meta.data, vars(all_of(resolutions))) %>% plyr::llply(levels)}) -> seurat

  # callModule(module=available_seurats.server, id='load_dataset')
  # seurat <- callModule(module=seurat_object.server, id='load_dataset')  # callModule(module=load_a_seurat.server, id='load_dataset')
  callModule(module=configuration_tab.server, id='configuration_tab', server_input=input, server_output=output, server_session=session)

  ## load the filter_seurat module
  # cell_filtering <- callModule(module=cell_filtering.server, id='seuratvis', seurat=seurat)

  ## load the servers for the analysis windows (menuItem or menuSubItem from the sidebar)
  callModule(module=highlight_feature_tab.server, id='highlight_feature_tab', server_input=input, server_output=output, server_session=session, seurat=seurat)
  callModule(module=highlight_feature_and_clusters_tab.server, id='highlight_feature_and_clusters_tab', server_input=input, server_output=output, server_session=session, seurat=seurat)
  callModule(module=gene_module_score_in_clusters_tab.server, id='gene_module_score_in_clusters_tab', server_input=input, server_output=output, server_session=session, seurat=seurat)
  callModule(module=provenance_tab.server, id='provenance_tab', server_input=input, server_output=output, server_session=session, seurat=seurat)

  observeEvent(input$simLoad, {
    seurat$object <- rnorm(1)
    seurat$object <- get('human_CS12', envir=globalenv())
    seurat$metadata <- seurat$object@meta.data
    seurat$gene_module_scores <- select_at(seurat$metadata, vars(starts_with('GeneModule-')))})

  # ###############################################################################################
  # reactions to tab selection
  ## react to opening tab with a filtered object loaded
  # observeEvent(input$sidebarmenu, {
  #   if(!is.null(seurat$object) & input$sidebarmenu=='cell_filtering-tab' && (!is.null(seurat$object@misc$cells_filtered) && seurat$object@misc$cells_filtered))
  #     sendSweetAlert(session=session, type='success', html=TRUE,
  #                    title='Notice', btn_labels='Great!',
  #                    text=tags$span('It looks like low-quality cells have already been removed from this Seurat object:', tags$h5(tags$code('@misc$cells_filtered == TRUE'))),
  #                    closeOnClickOutside=TRUE, showCloseButton=FALSE)})

  # ###############################################################################################
  # any code to exectue when the session ends
  session$onSessionEnded(function() {
    message('### session ended')})
}

