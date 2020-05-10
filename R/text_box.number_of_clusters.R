#' Display number of clusters in a dataset
#' 
#' Show the total number of clusters in a map
#' 
#' @param id unique name of the element
#' @param width integer value for width (rows sum to 12)
#' @param input,output,session used internally
#' 
#' @examples
#' 
#' \dontrun{
#' ui <- fluidPage(number_of_clusters_text_box.ui(id='page_name'))
#' server <- function(input, output, session) callModule(number_of_clusters_text_box.server, id='page_name')}
#' 
#' @rdname number_of_clusters_text_box
#' 
number_of_clusters_text_box.ui <- function(id, width=12) {
  module <- 'number_of_clusters'

  # make unique id for this object
  ns <- NS(namespace=id, id=module)
  module_env <- str_c(ns, 'env', sep='.')

  # create an environment in the seuratvis namespace
  e <- new.env()
  e$id <- id
  assign(x=module_env, val=e, envir=parent.frame(n=1))
  
  # return ui element(s)
  valueBoxOutput(outputId=ns, width=width)
}

#' @rdname number_of_clusters_text_box
#' 
number_of_clusters_text_box.server <- function(input, output, session) {
  # get environments containing variables to run/configure this object
  collect_environments(module='number_of_clusters') # provides `seuratvis_env` and `module_env`

  # make the text box
  renderValueBox(expr={
    # get the box subtitle
    switch(module_env$id,
           'Cell clusters') -> subtitle

    # create output object
    list(value={seuratvis_env$seurat_object.reactions$selected_clusters_per_resolution %>% comma()},
         subtitle=subtitle,
         icon=icon('first-order')) %>%
      modifyList(x=seuratvis:::text_box_defaults()) %>%
      do.call(what=valueBox)}) -> output$number_of_clusters
}
