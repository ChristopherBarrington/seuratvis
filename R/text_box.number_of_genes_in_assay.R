#' Display number of genes reported in a dataset assay
#' 
#' Show the total number of genes (rows) in the matrix of the active assay
#' 
#' @param id unique name of the element
#' @param width integer value for width (rows sum to 12)
#' @param input,oputput,session used internally
#' @param subtitle type of subtitle to display; one of: \code{default}
#' 
#' @examples
#' 
#' \dontrun{
#' ui <- fluidPage(number_of_genes_in_assay_text_box.ui(id='page_name'))
#' server <- function(input, output, session) callModule(number_of_genes_in_assay_text_box.server, id='page_name')}
#' 
#' @rdname number_of_genes_in_assay_text_box
#' 
number_of_genes_in_assay_text_box.ui <- function(id, width=12, subtitle='default') {
  module <- 'number_of_genes_in_assay'

  # make unique id for this object
  ns <- NS(namespace=id, id=module)
  module_env <- str_c(ns, 'env', sep='.')

  # create an environment in the seuratvis namespace
  e <- new.env()
  e$subtitle <- subtitle
  assign(x=module_env, val=e, envir=parent.frame(n=1))
  
  # return ui element(s)
  valueBoxOutput(outputId=ns, width=width)
}

#' @rdname number_of_genes_in_assay_text_box
#' 
number_of_genes_in_assay_text_box.server <- function(input, output, session) {
  # get environments containing variables to run/configure this object
  collect_environments(module='number_of_genes_in_assay') # provides `seuratvis_env` and `module_env`

  # make the text box
  renderValueBox(expr={
    switch(module_env$subtitle,
           default='Unique genes in assay',
           'Check the `subtitle` argument of your UI function!') -> subtitle

    list(value={nrow(seuratvis_env$seurat_object.reactions$seurat) %>% comma()},
         subtitle=subtitle,
         icon=icon('galactic-senate')) %>%
      modifyList(x=seuratvis:::text_box_defaults()) %>%
      do.call(what=valueBox)}) -> output$number_of_genes_in_assay
}
