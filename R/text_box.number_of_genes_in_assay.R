#' Display number of genes reported in a dataset assay
#' 
#' Show the total number of genes (rows) in the matrix of the active assay
#' 
#' @param id unique name of the element
#' @param width integer value for width (rows sum to 12)
#' @param input,output,session used internally
#' 
#' @examples
#' 
#' \dontrun{
#' ui <- fluidPage(number_of_genes_in_assay_text_box.ui(id='page_name'))
#' server <- function(input, output, session) callModule(number_of_genes_in_assay_text_box.server, id='page_name')}
#' 
#' @rdname number_of_genes_in_assay_text_box
#' 
number_of_genes_in_assay_text_box.ui <- function(id, width=12) {
  sprintf(fmt='### %s-number_of_genes_in_assay_text_box.ui', id) %>% message()

  module <- 'number_of_genes_in_assay_text_box'

  # make unique id for this object
  ns <- NS(namespace=id, id=module)

  # create an environment in the seuratvis namespace
  e <- new.env()
  e$id <- id
  assign(x=ns, val=e, envir=module_environments)

  # record the server(s) to call
  get0(env=module_servers_to_call, x=id) %>% append(sprintf(fmt='%s.server', module)) %>% assign(env=module_servers_to_call, x=id)

  # return ui element(s)
  valueBoxOutput(outputId=ns, width=width)
}

#' @rdname number_of_genes_in_assay_text_box
#' 
number_of_genes_in_assay_text_box.server <- function(input, output, session) {
  session$ns('') %>% sprintf(fmt='### %snumber_of_genes_in_assay_text_box.server') %>% message()

  # get environments containing variables to run/configure this object
  collect_environments(id=parent.frame()$id, module='number_of_genes_in_assay_text_box') # provides `seuratvis_env`, `server_env` and `module_env`

  # make the text box
  renderValueBox(expr={
    # get the box subtitle
    switch(module_env$id,
           'Unique features in assay') -> subtitle

    # create output object
    list(value={seurat_object.reactions$reference_metrics$n_features %>% comma()},
         subtitle=subtitle,
         icon=icon('galactic-senate')) %>%
      modifyList(x=seuratvis:::text_box_defaults()) %>%
      do.call(what=valueBox)}) -> output$number_of_genes_in_assay_text_box
}
