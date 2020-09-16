#' 
#' 
show_filtering_parameters.ui <- function(id, label='Cell filtering parameters', include_copy_buttons=TRUE) {
  # make ui elements
  ## a text output box
  ace_editor.ui(id=NS(id, 'ace_verbatim_text_output'), showLineNumbers=FALSE, theme='tomorrow_night_eighties') -> ace_text_output_box

  ## copy buttons are required
  div(tags$h6('Copy to clipboard'),
      uiOutput(outputId=NS(id, 'plain.copybutton'), inline=TRUE),
      uiOutput(outputId=NS(id, 'csv.copybutton'), inline=TRUE),
      uiOutput(outputId=NS(id, 'r.copybutton'), inline=TRUE)) -> copy_buttons

  # return ui element(s)
  tagList(tags$label(label), ace_text_output_box, if(include_copy_buttons) copy_buttons)
}

#' 
#' 
show_filtering_parameters.server <- function(input, output, session, seurat, cell_filtering, filters=list()) {
  filtering_arguments <- reactiveValues()

  display_text <- reactiveVal()
  callModule(module=ace_editor.server, id='ace_verbatim_text_output', display_text=display_text)

  # react to changes in the filters
  observe({
    req(seurat$project)

    # if there are no filters provided
    if(length(filters)==0)
      return(NULL)

    # if any of the filters are not uet initialised
    if(lapply(filters, reactiveValuesToList) %>% sapply(length) %>% equals(0) %>% any())
      return(NULL)

    # get the values in the list of reactives
    lapply(filters, reactiveValuesToList)  %>%
      lapply(function(x) x %>% extract(str_detect(string=names(x), pattern='variable|min|max|in_set'))) %>% # pick out the elements we can use
      plyr::ldply(as.data.frame) %>%
      gather(key=logic, value=value, -variable) %>%
      drop_na() -> filters_df

    # prepare the condition filter(s)
    all_subset_conditions <- NULL
    if(nrow(filters_df)>0)
      filters_df %>%
        mutate(logic=factor(logic, levels=c('min', 'max', 'in_set')),
               logic=fct_recode(logic, `>=`='min', `<=`='max', ` %in% `='in_set'),
               value=str_trim((value))) %>%
        arrange(variable, value) %>%
        apply(1, str_c, collapse='') -> all_subset_conditions

    # combine all output lines
    list(project_line={seurat$project %>% sprintf(fmt='# %s')},
         n_cells_line={cell_filtering$n_cells %>% comma() %>% sprintf(fmt='# n_cells=%s')},
         n_umi_line={cell_filtering$n_umi %>% comma() %>% sprintf(fmt='# n_umi=%s')},
         filters_line={all_subset_conditions %>% str_c(collapse=' &\n')}) %>%
      unlist() %>%
      str_c(collapse='\n') -> output_text

    # update the ui with filtering parameters
    display_text('refresh')
    display_text(output_text)

    # update the reactive
    filtering_arguments$output_text <- output_text
    filtering_arguments$all_subset_conditions <- all_subset_conditions})

  # prepare copy to clipboard buttons
  ## copy text as is
  renderUI(expr={
    filtering_arguments$output_text %>%
      str_c('\n') %>%
      rclipButton(inputId='rclipButton.plain.in', label='', icon('clipboard-check'))}) -> output$plain.copybutton

  ## comma separate elements
  renderUI(expr={
    c(seurat$project, filtering_arguments$all_subset_conditions) %>%
      str_c(collapse=',') %>%
      str_c('\n') %>%
      rclipButton(inputId='rclipButton.csv.in', label='', icon('file-excel'))}) -> output$csv.copybutton

  ## copy only the conditional
  renderUI(expr={
    filtering_arguments$all_subset_conditions %>%
      str_c(collapse=' & ') %>%
      str_c('\n') %>%
      rclipButton(inputId='rclipButton.r.in', label='', icon('r-project'))}) -> output$r.copybutton
}
