#' Display project name
#' 
#' Create a formatted project name text box
#' 
#' @param id unique name of the element
#' @param width integer value for width (rows sum to 12)
#' @param input,output,session used internally
#' 
#' @examples
#' 
#' \dontrun{
#' ui <- fluidPage(project_name_text_box.ui(id='page_name'))
#' server <- function(input, output, session) callModule(project_name_text_box.server, id='page_name')}
#' 
#' @rdname project_name_text_box
#' 
project_name_text_box.ui <- function(id, width=12) {
  sprintf(fmt='### %s-project_name_text_box.ui', id) %>% message()

  module <- 'project_name_text_box'

  # record the server(s) to call
  get0(env=module_servers_to_call, x=id) %>% append(sprintf(fmt='%s.server', module)) %>% assign(env=module_servers_to_call, x=id)

  # return ui element(s)
  valueBoxOutput(outputId=NS(namespace=id, id=module), width=width)
}

#' @rdname project_name_text_box
#' 
project_name_text_box.server <- function(input, output, session, seurat, ...) {
  session$ns('') %>% sprintf(fmt='### %sproject_name_text_box.server') %>% message()

  renderValueBox(expr={
    # send a message
    session$ns('') %>% sprintf(fmt='### %sproject_name_text_box.server-renderValueBox') %>% message()

    list(value=seurat$formatted_project,
         subtitle='Loaded Seurat object',
         icon=icon('certificate')) %>%
    modifyList(x=text_box_defaults()) %>%
    do.call(what=valueBox)}) -> output$project_name_text_box
}
