#' Display median number of reads per cell
#' 
#' Show the median number of reads in cells of the dataset
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
#' ui <- fluidPage(number_of_reads_per_cell_text_box.ui(id='page_name'))
#' server <- function(input, output, session) callModule(number_of_reads_per_cell_text_box.server, id='page_name')}
#' 
#' @rdname median_reads_per_cell_text_box
#' 
number_of_reads_per_cell_text_box.ui <- function(id, width=12, subtitle='default', f=median, f_name=deparse(substitute(f))) {
  module <- 'number_of_reads_per_cell'

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

#' @rdname median_reads_per_cell_text_box
#' 
number_of_reads_per_cell_text_box.server <- function(input, output, session) {
  # get environments containing variables to run/configure this object
  collect_environments(module='number_of_reads_per_cell') # provides `seuratvis_env` and `module_env`

  # make the text box
  renderValueBox(expr={
    switch(module_env$subtitle,
           default=sprintf('%s reads per cell', str_to_title(module_env$summary_function_name)),
           'Check the `subtitle` argument of your UI function!') -> subtitle

    list(value={seuratvis_env$seurat_object.reactions$seurat$nCount_RNA %>% module_env$summary_function() %>% round(digits=1) %>% comma()},
         subtitle=subtitle,
         icon=icon('frog')) %>%
      modifyList(x=seuratvis:::text_box_defaults()) %>%
      do.call(what=valueBox)}) -> output$number_of_reads_per_cell
}
