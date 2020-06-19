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
show_filtering_parameters.server <- function(input, output, session, seurat, cell_filtering, ...) {
  session$ns('') %>% sprintf(fmt='### %sshow_filtering_parameters.server') %>% message()

  # get environments containing variables to run/configure this object
  collect_environments(id=parent.frame()$id, module='show_filtering_parameters') # provides `seuratvis_env`, `server_env` and `module_env`
  session_server <- get(x='session', env=server_env)

  format_subset_conditional <- function(x, fmt) ifelse(is.na(x), NA, sprintf(fmt=fmt, x))
  group_format_subset_conditional <- function(x) x %>% na.omit() %>% paste(collapse=' & ')
  filtering_arguments <- list()

  # observeEvent(eventExpr=reactiveValuesToList(filtering_parameters.reactions), handlerExpr={
  observeEvent(eventExpr=c(seurat$done, cell_filtering$done), handlerExpr={
    # make sure these elements are defined
    req(seurat$n_features_variable)
    req(seurat$n_umi_variable)
    req(seurat$proportion_mt_variable)

    # send a message
    sprintf('### %sshow_filtering_parameters.server-observeEvent-cell_filtering$done', session$ns('')) %>% message()

    # save formatted filters
    c(format_subset_conditional(x=cell_filtering$total_umi_per_cell_min, fmt='X>=%d'),
      format_subset_conditional(x=cell_filtering$total_umi_per_cell_max, fmt='X<=%d')) %>%
      str_replace_all(pattern='X', replacement=seurat$n_umi_variable) %>%
      group_format_subset_conditional() -> umi_filter

    c(format_subset_conditional(x=cell_filtering$features_per_cell_min, fmt='X>=%d'),
      format_subset_conditional(x=cell_filtering$features_per_cell_max, fmt='X<=%d')) %>%
      str_replace_all(pattern='X', replacement=seurat$n_features_variable) %>%
      group_format_subset_conditional() -> features_filter

    format_subset_conditional(x=cell_filtering$max_percent_mitochondria, fmt='X<=%s') %>%
      str_replace_all(pattern='X', replacement=seurat$proportion_mt_variable) -> mt_filter

    all_subset_conditions <- list(umi_filter, features_filter, mt_filter)
    filtering_arguments$all_subset_conditions <- all_subset_conditions

    # combine all output lines
    list(project_line={seurat$project %>% sprintf(fmt='# %s')},
         n_cells_line={cell_filtering$n_cells %>% comma() %>% sprintf(fmt='# n_cells=%s')},
         n_umi_line={cell_filtering$n_umi %>% comma() %>% sprintf(fmt='# n_umi=%s')},
         filters_line={all_subset_conditions %>% str_c(collapse=' &\n')}) %>%
      str_c(collapse='\n') -> output_text

    # update the ui with filtering parameters
    renderText({output_text}) -> output$verbatim_text_output
    filtering_arguments$output_text <- output_text})

  # prepare copy to clipboard buttons
  ## copy text as is
  renderUI(expr={
    rclipButton(inputId='rclipButton.plain.in', label='', icon('clipboard-check'),
                clipText=filtering_arguments$output_text)}) -> output$plain.copybutton
  ## comma separate elements
  renderUI(expr={
    c(seurat$project,
      paste(filtering_arguments$all_subset_conditions, collapse=' & ')) %>%
      paste(collapse=',') %>%
      rclipButton(inputId='rclipButton.csv.in', label='', icon('file-excel'))}) -> output$csv.copybutton

  ## copy only the conditional
  renderUI(expr={
    paste(filtering_arguments$all_subset_conditions, collapse=' & ') %>%
      rclipButton(inputId='rclipButton.r.in', label='', icon('r-project'))}) -> output$r.copybutton
}
