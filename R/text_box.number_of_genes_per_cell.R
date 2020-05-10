#' Display number of genes per cell
#' 
#' Show the average number of genes in cells of the dataset
#' 
#' @param id unique name of the element
#' @param width integer value for width (rows sum to 12)
#' @param input,oputput,session used internally
#' @param subtitle type of subtitle to display; one of: \code{default}
#' @param f,f_name summary function taking one vector and returning a value and a name used to make the subtitle
#' 
#' @examples
#' 
#' \dontrun{
#' ui <- fluidPage(number_of_genes_per_cell_text_box.ui(id='page_name'))
#' server <- function(input, output, session) callModule(number_of_genes_per_cell_text_box.server, id='page_name')}
#' 
#' @rdname number_of_genes_per_cell_text_box
#' 
number_of_genes_per_cell_text_box.ui <- function(id, width=12, subtitle='default', f=median, f_name=deparse(substitute(f))) {
  module <- 'number_of_genes_per_cell'

  # make unique id for this object
  ns <- NS(namespace=id, id=module)
  module_env <- str_c(ns, 'env', sep='.')

  # create an environment in the seuratvis namespace
  e <- new.env()
  e$subtitle <- subtitle
  e$summary_function <- f
  e$summary_function_name <- f_name
  assign(x=module_env, val=e, envir=parent.frame(n=1))

  # return ui element(s)
  valueBoxOutput(outputId=ns, width=width)
}

#' @rdname number_of_genes_per_cell_text_box
#' 
number_of_genes_per_cell_text_box.server <- function(input, output, session) {
  # get environments containing variables to run/configure this object
  collect_environments(module='number_of_genes_per_cell') # provides `seuratvis_env` and `module_env`

  # make the text box
  renderValueBox(expr={
    switch(module_env$id,
           cell_filtering=sprintf(fmt='%s genes per cell (%s)', str_to_title(module_env$summary_function_name), comma(seuratvis_env$cell_filtering_data.reactions$median_genes_per_cell-seuratvis_env$seurat_object.reactions$reference_metrics$median_genes_per_cell) %>% ifelse(str_detect(., '^-'), ., str_c('+', .))),
           sprintf('%s genes per cell', str_to_title(module_env$summary_function_name))) -> subtitle

    list(value={seuratvis_env$cell_filtering_data.reactions$filtered_cell_set$nFeature_RNA %>% module_env$summary_function() %>% round(digits=1) %>% comma()},
         subtitle=subtitle,
         icon=icon('crow')) %>%
      modifyList(x=seuratvis:::text_box_defaults()) %>%
      do.call(what=valueBox)}) -> output$number_of_genes_per_cell
}
