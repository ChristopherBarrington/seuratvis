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
  sprintf(fmt='### %s-show_filtering_parameters.ui', id) %>% message()

  module <- 'show_filtering_parameters'

  # make unique id for this object
  ns <- NS(namespace=id)
  module_ns <- ns(id=module)

  # create an environment in the seuratvis namespace
  e <- new.env()
  e$id <- id
  assign(x=module_ns, val=e, envir=module_environments)

  # module_environments$show_filtering_parameters$ns %<>% c(module_ns)
  # module_environments$show_filtering_parameters$id %<>% c(id)

  # record the server(s) to call
  get0(env=module_servers_to_call, x=id) %>% append(sprintf(fmt='%s.server', module)) %>% assign(env=module_servers_to_call, x=id)

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
  session$ns('') %>% sprintf(fmt='### %sshow_filtering_parameters.server') %>% message()

  # get environments containing variables to run/configure this object
  collect_environments(id=parent.frame()$id, module='show_filtering_parameters') # provides `seuratvis_env`, `server_env` and `module_env`
  session_server <- get(x='session', env=server_env)

  format_subset_conditional <- function(x, fmt) ifelse(is.na(x), NA, sprintf(fmt=fmt, x))
  group_format_subset_conditional <- function(x) x %>% na.omit() %>% paste(collapse=' & ')

  observeEvent(eventExpr=reactiveValuesToList(filtering_parameters.reactions), handlerExpr={
    # make sure required reactives are available
    req(seurat_object.reactions$project)
    req(seurat_configuration.reactions$n_features_variable)
    req(seurat_configuration.reactions$n_umi_variable)
    req(seurat_configuration.reactions$proportion_mt_variable)
    req(filtered_cells.reactions$n_cells)

    message('### show_filtering_parameters.server-observeEvent-reactiveValuesToList(filtering_parameters.reactions)')

    # create variables for shorthand
    thresholds <- reactiveValuesToList(filtering_parameters.reactions)
    filtered_cells <- reactiveValuesToList(filtered_cells.reactions)

    # save formatted filters
    c(format_subset_conditional(x=thresholds$total_umi_per_cell_min, fmt='X>=%d'),
      format_subset_conditional(x=thresholds$total_umi_per_cell_max, fmt='X<=%d')) %>%
      str_replace_all(pattern='X', replacement=seurat_configuration.reactions$n_umi_variable) %>%
      group_format_subset_conditional() -> umi_filter

    c(format_subset_conditional(x=thresholds$features_per_cell_min, fmt='X>=%d'),
      format_subset_conditional(x=thresholds$features_per_cell_max, fmt='X<=%d')) %>%
      str_replace_all(pattern='X', replacement=seurat_configuration.reactions$n_features_variable) %>%
      group_format_subset_conditional() -> features_filter

    format_subset_conditional(x=thresholds$max_percent_mitochondria, fmt='X<=%s') %>%
      str_replace_all(pattern='X', replacement=seurat_configuration.reactions$proportion_mt_variable) -> mt_filter

    all_subset_conditions <- list(umi_filter, features_filter, mt_filter)
    filtering_arguments.reactions$all_subset_conditions <- all_subset_conditions

    # combine all output lines
    list(project_line={seurat_object.reactions$project %>% sprintf(fmt='# %s')},
         n_cells_line={filtered_cells$n_cells %>% comma() %>% sprintf(fmt='# n_cells=%s')},
         n_umi_line={filtered_cells$n_umi %>% comma() %>% sprintf(fmt='# n_umi=%s')},
         filters_line={all_subset_conditions %>% str_c(collapse=' &\n')}) %>%
      str_c(collapse='\n') -> output_text

    # update the ui with filtering parameters
    renderText({output_text}) -> output$verbatim_text_output
    filtering_arguments.reactions$output_text <- output_text})

  # prepare copy to clipboard buttons
  ## copy text as is
  renderUI(expr={
    rclipButton(inputId='rclipButton.plain.in', label='', icon('clipboard-check'),
                clipText=filtering_arguments.reactions$output_text)}) -> output$plain.copybutton
  ## comma separate elements
  renderUI(expr={
    c(filtering_parameters.reactions$project,
      paste(filtering_arguments.reactions$all_subset_conditions, collapse=' & ')) %>%
      paste(collapse=',') %>%
      rclipButton(inputId='rclipButton.csv.in', label='', icon('file-excel'))}) -> output$csv.copybutton

  ## copy only the conditional
  renderUI(expr={
    paste(filtering_arguments.reactions$all_subset_conditions, collapse=' & ') %>%
      rclipButton(inputId='rclipButton.r.in', label='', icon('r-project'))}) -> output$r.copybutton
}
