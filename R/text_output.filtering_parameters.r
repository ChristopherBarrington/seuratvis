#' 
#' 
show_filtering_parameters.ui <- function(id, label='Cell filtering parameters', include_copy_buttons=TRUE) {
  # make ui elements
  ## a text output box
  verbatimTextOutput(outputId=NS(id, 'verbatim_text_output'), placeholder=TRUE) -> text_output_box

  ## copy buttons are required
  div(tags$h6('Copy to clipboard'),
      uiOutput(outputId=NS(id, 'plain.copybutton'), inline=TRUE),
      uiOutput(outputId=NS(id, 'csv.copybutton'), inline=TRUE),
      uiOutput(outputId=NS(id, 'r.copybutton'), inline=TRUE)) -> copy_buttons

  # return ui element(s)
  tagList(tags$label(label), text_output_box, if(include_copy_buttons) copy_buttons)
}

#' 
#' 
show_filtering_parameters.server <- function(input, output, session, seurat, cell_filtering) {
  format_subset_conditional <- function(x, fmt) ifelse(is.na(x), NA, sprintf(fmt=fmt, x))
  group_format_subset_conditional <- function(x, sep=' & ') x %>% na.omit() %>% paste(collapse=sep)
  filtering_arguments <- reactiveValues()

  # observeEvent(eventExpr=reactiveValuesToList(filtering_parameters.reactions), handlerExpr={
  observeEvent(eventExpr=cell_filtering$done, handlerExpr={
    # make sure these elements are defined
    req(seurat$n_features_variable)
    req(seurat$n_umi_variable)
    req(seurat$proportion_mt_variable)

    # save formatted filters
    c(format_subset_conditional(x=cell_filtering$n_umi_min, fmt='X>=%d'),
      format_subset_conditional(x=cell_filtering$n_umi_max, fmt='X<=%d')) %>%
      str_replace_all(pattern='X', replacement=seurat$n_umi_variable) -> umi_filter

    c(format_subset_conditional(x=cell_filtering$n_features_min, fmt='X>=%d'),
      format_subset_conditional(x=cell_filtering$n_features_max, fmt='X<=%d')) %>%
      str_replace_all(pattern='X', replacement=seurat$n_features_variable) -> features_filter

    format_subset_conditional(x=cell_filtering$proportion_mt_max, fmt='X<=%s') %>%
      str_replace_all(pattern='X', replacement=seurat$proportion_mt_variable) -> mt_filter

    all_subset_conditions <- list(umi_filter=umi_filter, features_filter=features_filter, mt_filter=mt_filter)
    filtering_arguments$all_subset_conditions <- all_subset_conditions

    # combine all output lines
    list(project_line={seurat$project %>% sprintf(fmt='# %s')},
         n_cells_line={cell_filtering$n_cells %>% comma() %>% sprintf(fmt='# n_cells=%s')},
         n_umi_line={cell_filtering$n_umi %>% comma() %>% sprintf(fmt='# n_umi=%s')},
         filters_line={all_subset_conditions %>% unlist() %>% str_c(collapse=' &\n')}) %>%
      str_c(collapse='\n') -> output_text

    # update the ui with filtering parameters
    renderText({output_text}) -> output$verbatim_text_output
    filtering_arguments$output_text <- output_text})

  # prepare copy to clipboard buttons
  ## copy text as is
  renderUI(expr={
    filtering_arguments$output_text %>%
      str_c('\n') %>%
      rclipButton(inputId='rclipButton.plain.in', label='', icon('clipboard-check'))}) -> output$plain.copybutton

  ## comma separate elements
  renderUI(expr={
    c(seurat$project,
      unlist(filtering_arguments$all_subset_conditions)) %>%
      str_c(collapse=',') %>%
      str_c('\n') %>%
      rclipButton(inputId='rclipButton.csv.in', label='', icon('file-excel'))}) -> output$csv.copybutton

  ## copy only the conditional
  renderUI(expr={
    filtering_arguments$all_subset_conditions %>%
      unlist() %>%
      str_c(collapse=' & ') %>%
      str_c('\n') %>%
      rclipButton(inputId='rclipButton.r.in', label='', icon('r-project'))}) -> output$r.copybutton
}
