#' Print specified filtering parameters
#' 
#' Displays the currently applied filters and can provide formatted copy buttons
#' 
#' @param id unique name of the element
#' @param label text label of the element
#' @param include_copy_buttons should buttons to copy formatted text be included?
#' 
#' @examples
#' 
#' \dontrun{
#' ui <- fluidPage(show_filtering_parameters.ui(id='page_name'))
#' server <- function(input, output, session) {
#'   callModule(show_filtering_parameters.server, id='page_name')
#' }}
#' 
#' #' @rdname total_umi_per_cell_filter
#' 
show_filtering_parameters.ui <- function(id, label='Cell filtering parameters', include_copy_buttons=TRUE) {
  message('### show_filtering_parameters.ui')

  module <- 'show_filtering_parameters'

  # make unique id for this object
  ns <- NS(namespace=id)
  module_ns <- ns(id=module)

  # create an environment in the seuratvis namespace
  e <- new.env()
  e$id <- id
  assign(x=module_ns, val=e, envir=module_environments)

  module_environments$show_filtering_parameters$ns %<>% c(module_ns)
  module_environments$show_filtering_parameters$id %<>% c(id)

  # make ui elements
  ## if a label switch is required, make one
  verbatimTextOutput(outputId=ns(id='verbatim_text_output'), placeholder=TRUE) -> text_output_box

  ## if copy buttons are required
  copy_buttons <- NULL
  if(include_copy_buttons)
    div(tags$h6('Copy to clipboard', style='display: inline;'),
        uiOutput(outputId=ns(id='plain.copybutton'), inline=TRUE),
        uiOutput(outputId=ns(id='csv.copybutton'), inline=TRUE),
        uiOutput(outputId=ns(id='r.copybutton'), inline=TRUE)) -> copy_buttons

  # return ui element(s)
  tagList(tags$label(label), text_output_box, copy_buttons)
}

#' React to filtering parameters
#' 
#' @details
#' Updates the filtering parameters reactive values.
#' 
#' @rdname show_filtering_parameters
#' 
show_filtering_parameters.server <- function(input, output, session) {
  message('### show_filtering_parameters.server')

  # get environments containing variables to run/configure this object
  collect_environments(id=parent.frame()$id, module='show_filtering_parameters') # provides `seuratvis_env`, `server_env` and `module_env`
  session_server <- get(x='session', env=server_env)

  format_subset_conditional <- function(x, fmt) ifelse(is.na(x), NA, sprintf(fmt=fmt, x))
  group_format_subset_conditional <- function(x) x %>% na.omit() %>% paste(collapse=' & ')

  observeEvent(eventExpr=reactiveValuesToList(filtering_parameters.reactions), handlerExpr={
    message('### show_filtering_parameters.server-observeEvent-seurat_object.reactions$seurat')

    # make sure seurat object is loaded
    req(seurat_object.reactions$seurat)

    # create variables for shorthand
    thresholds <- reactiveValuesToList(filtering_parameters.reactions)

    # save formatted filters
    c(format_subset_conditional(x=thresholds$total_umi_per_cell_min, fmt='nCount_RNA>=%d'),
      format_subset_conditional(x=thresholds$total_umi_per_cell_max, fmt='nCount_RNA<=%d')) %>%
      group_format_subset_conditional() -> umi_filter

    c(format_subset_conditional(x=thresholds$features_per_cell_min, fmt='nFeature_RNA>=%d'),
      format_subset_conditional(x=thresholds$features_per_cell_max, fmt='nFeature_RNA<=%d')) %>%
      group_format_subset_conditional() -> features_filter

    'percent_mt' %>%
      sprintf(fmt='%s<=%%s') %>%
      format_subset_conditional(x=thresholds$max_percent_mitochondria) -> mt_filter

    all_subset_conditions <- list(umi_filter, features_filter, mt_filter)
    filtering_parameters.reactions$all_subset_conditions <- all_subset_conditions

    # combine all output lines
    list(project_line={thresholds$project %>% sprintf(fmt='# %s')},
         n_cells_line={thresholds$n_cells %>% comma() %>% sprintf(fmt='# n_cells=%s')},
         n_umi_line={thresholds$n_umi %>% comma() %>% sprintf(fmt='# n_umi=%s')},
         filters_line={all_subset_conditions %>% str_c(collapse=' &\n')}) %>%
      str_c(collapse='\n') -> output_text

    # update the ui with filtering parameters
    renderText({output_text}) -> output$verbatim_text_output
    filtering_parameters.reactions$output_text <- output_text})

  # prepare copy to clipboard buttons
  ## copy text as is
  renderUI(expr={
    rclipButton(inputId='rclipButton.plain.in', label='', icon('clipboard-check'),
                clipText=filtering_parameters.reactions$output_text)}) -> output$plain.copybutton
  ## comma separate elements
  renderUI(expr={
    c(filtering_parameters.reactions$project,
      paste(filtering_parameters.reactions$all_subset_conditions, collapse=' & ')) %>%
      paste(collapse=',') %>%
      rclipButton(inputId='rclipButton.csv.in', label='', icon('file-excel'))}) -> output$csv.copybutton

  ## copy only the conditional
  renderUI(expr={
    paste(filtering_parameters.reactions$all_subset_conditions, collapse=' & ') %>%
      rclipButton(inputId='rclipButton.r.in', label='', icon('r-project'))}) -> output$r.copybutton

  # initialise filtering reactive when Seurat object is loaded
  observeEvent(eventExpr=seurat_object.reactions$seurat, handlerExpr={
    message('### show_filtering_parameters.server-observeEvent-seurat_object.reactions$seurat')

    # create variables for shorthand
    seurat <- seurat_object.reactions$seurat

    list(project=Project(seurat),
         n_cells=nrow(seurat@meta.data),
         n_umi=sum(seurat@meta.data$nCount_RNA),
         
         total_umi_per_cell_min=min(seurat@meta.data$nCount_RNA),
         total_umi_per_cell_max=max(seurat@meta.data$nCount_RNA),
         
         features_per_cell_min=min(seurat@meta.data$nFeature_RNA),
         features_per_cell_max=max(seurat@meta.data$nFeature_RNA),
         
         max_percent_mitochondria=round(max(seurat@meta.data$percent_mt)+0.05, digits=1)) -> initial_values

    for(i in names(initial_values))
      filtering_parameters.reactions[[i]] <- initial_values[[i]]})
}
