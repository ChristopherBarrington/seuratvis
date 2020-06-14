#' Display number of cells in a dataset
#' 
#' Show the total number of cells
#' 
#' @param id unique name of the element
#' @param width integer value for width (rows sum to 12)
#' @param input,output,session used internally
#' 
#' @examples
#' 
#' \dontrun{
#' ui <- fluidPage(number_of_cells_text_box.ui(id='page_name'))
#' server <- function(input, output, session) callModule(number_of_cells_text_box.server, id='page_name')}
#' 
#' @rdname number_of_cells_text_box
#' 
number_of_cells_text_box.ui <- function(id, width=12) {
  sprintf(fmt='### %s-number_of_cells_text_box.ui', id) %>% message()

  module <- 'number_of_cells_text_box'

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

#' @import scales
#' 
#' @rdname number_of_cells_text_box
#' 
number_of_cells_text_box.server <- function(input, output, session, seurat, cell_filtering, ...) {
  session$ns('') %>% sprintf(fmt='### %snumber_of_cells_text_box.server') %>% message()

  # get environments containing variables to run/configure this object
  collect_environments(id=parent.frame()$id, module='number_of_cells_text_box') # provides `seuratvis_env`, `server_env` and `module_env`

  # make the text box
  renderValueBox(expr={
    # send a message
    session$ns('') %>% sprintf(fmt='### %snumber_of_cells_text_box.server-renderValueBox') %>% message()

    # create variables for shorthand
    n_reference <- seurat$n_cells
    n_filtered <- cell_filtering$n_cells

    # get the box subtitle
    #! TODO: use a type argument to ui function to make separate servers and avoid using module_env?
    switch(module_env$id,
           cell_filtering=sprintf(fmt='Cells remaining (%.1f%%)', n_filtered/n_reference*100),
           'Cells in map') -> subtitle

    # create output object
    list(value={n_filtered %>% comma()},
         subtitle=subtitle,
         icon=icon('galactic-republic')) %>%
      modifyList(x=seuratvis:::text_box_defaults()) %>%
      do.call(what=valueBox)}) -> output$number_of_cells_text_box
}
