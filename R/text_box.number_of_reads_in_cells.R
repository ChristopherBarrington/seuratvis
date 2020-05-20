#' Display number of reads in a dataset
#' 
#' Show the total number of reads in cells
#' 
#' @param id unique name of the element
#' @param width integer value for width (rows sum to 12)
#' @param input,output,session used internally
#' 
#' @examples
#' 
#' \dontrun{
#' ui <- fluidPage(number_of_reads_text_box.ui(id='page_name'))
#' server <- function(input, output, session) callModule(number_of_reads_text_box.server, id='page_name')}
#' 
#' @rdname number_of_reads_text_box
#' 
number_of_reads_text_box.ui <- function(id, width=12) {
  module <- 'number_of_reads'

  # make unique id for this object
  ns <- NS(namespace=id, id=module)

  # create an environment in the seuratvis namespace
  e <- new.env()
  e$id <- id
  assign(x=ns, val=e, envir=module_environments)
  
  # return ui element(s)
  valueBoxOutput(outputId=ns, width=width)
}

#' @rdname number_of_reads_text_box
#' 
number_of_reads_text_box.server <- function(input, output, session) {
  # get environemtns containing variables to run/configure this object
  collect_environments(id=parent.frame()$id, module='number_of_reads') # provides `seuratvis_env`, `server_env` and `module_env`

  # make the text box
  renderValueBox(expr={
    server_env$react_to_cell_filtering() # not sure about this line, is it necessary?

    # get the box subtitle
    switch(module_env$id,
           cell_filtering=sprintf(fmt='Total reads remaining (%.1f%%)', server_env$cell_filtering_data.reactions$total_reads/seurat_object.reactions$reference_metrics$total_reads*100), # neaten this up
           'Total reads in cells') -> subtitle

    # create output object
    list(value={server_env$cell_filtering_data.reactions$total_reads %>% comma()}, # not the best, should not depend on this object, may not be filtered
         subtitle=subtitle,
         icon=icon('old-republic')) %>%
      modifyList(x=seuratvis:::text_box_defaults()) %>%
      do.call(what=valueBox)}) -> output$number_of_reads
}
