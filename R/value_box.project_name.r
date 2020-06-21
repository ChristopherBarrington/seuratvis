#' 
#' 
project_name_text_box.ui <- function(id, width=12)
  valueBoxOutput(outputId=NS(id, 'box'), width=width)

#'
#'  
project_name_text_box.server <- function(input, output, session, seurat) {
  renderValueBox(expr={
    list(value=seurat$formatted_project_name,
         subtitle='Loaded Seurat object',
         icon=icon('certificate')) %>%
    modifyList(x=list(color='purple')) %>%
    do.call(what=valueBox)}) -> output$box
}
