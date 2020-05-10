#' Display number of reads in a dataset
#' 
#' Show the total number of reads in cells
#' 
#' @param id unique name of the element
#' @param width integer value for width (rows sum to 12)
#' @param input,oputput,session used internally
#' @param subtitle type of subtitle to display; one of: \code{default}, \code{proportion}
#' 
#' @examples
#' 
#' \dontrun{
#' ui <- fluidPage(number_of_reads_text_box.ui(id='page_name'))
#' server <- function(input, output, session) callModule(number_of_reads_text_box.server(id='page_name'))}
#' 
#' @rdname number_of_reads_text_box
#' 
number_of_reads_text_box.ui <- function(id, width=12, subtitle='default') {
  module <- 'number_of_reads'

  # make unique id for this object
  ns <- NS(namespace=id, id=module)
  module_env <- str_c(ns, 'env', sep='.')

  # create an environment in the seuratvis namespace
  e <- new.env()
  e$subtitle <- subtitle
  assign(x=module_env, val=e, envir=parent.frame(n=1)) # don't like using parent.frame, can use anything better?
  
  # return ui element(s)
  valueBoxOutput(outputId=ns, width=width)
}

#' @rdname number_of_reads_text_box
#' 
number_of_reads_text_box.server <- function(input, output, session) {
  # get environemtns containing variables to run/configure this object
  collect_environments(module='number_of_reads') # provides `seuratvis_env` and `module_env`

  # make the text box
  renderValueBox(expr={
    seuratvis_env$react_to_cell_filtering() # not sure about this line, is it necessary?
    switch(module_env$subtitle,
           proportion=sprintf(fmt='Total reads remaining (%.1f%%)', seuratvis_env$cell_filtering_data.reactions$total_reads/seuratvis_env$seurat_object.reactions$reference_metrics$total_reads*100), # neaten this up
           default='Total reads in cells',
           'Check the `subtitle` argument of your UI function!') -> subtitle

    list(value=scales::comma(seuratvis_env$cell_filtering_data.reactions$total_reads), # not the best, should not depend on this object, may not be filtered
         subtitle=subtitle,
         icon=icon('old-republic')) %>%
      modifyList(x=seuratvis:::text_box_defaults()) %>%
      do.call(what=valueBox)}) -> output$number_of_reads
}
