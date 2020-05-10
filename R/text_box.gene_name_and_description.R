#' Display gene information
#' 
#' Displays gene of interest and biomaRt description
#' 
#' @param id unique name of the element
#' @param width integer value for width (rows sum to 12)
#' @param input,oputput,session used internally
#' 
#' @examples
#' 
#' \dontrun{
#' ui <- fluidPage(gene_name_and_description_text_box.ui(id='page_name'))
#' server <- function(input, output, session) callModule(gene_name_and_description_text_box.server, id='page_name')}
#' 
#' @rdname gene_name_and_description_text_box
#' 
gene_name_and_description_text_box.ui <- function(id, width=12)
  valueBoxOutput(outputId=NS(namespace=id, id='gene_name_and_description'), width=width)

#' @rdname gene_name_and_description_text_box
#' 
gene_name_and_description_text_box.server <- function(input, output, session) {
  renderValueBox(env=parent.frame(n=2), quoted=FALSE, expr={
    Misc(seurat_object.reactions$seurat, slot='mart') %>%
      get_formatted_gene_description(external_gene_name=input$gene_of_interest.dd) -> subtitle

    list(value=input$gene_of_interest.dd,
         subtitle=subtitle,
         icon=icon('jedi-order')) %>%
      modifyList(x=text_box_defaults()) %>%
      do.call(what=valueBox)}) -> output$gene_name_and_description
}

#' Get description of a gene
#' 
#' Use biomaRt to get formatted description of a gene from gene name
#' 
#' @param external_gene_name name of gene to query
#' @param mart a \code{biomaRt} object
#' 
#' @importFrom biomaRt getBM
#' 
#' @return a character string of the \code{attribute} from the \code{mart} for the \code{external_gene_name}.
#' 
get_formatted_gene_description <- function(external_gene_name, mart) {
  tryCatch(expr=getBM(mart=mart, attributes='description', filter='external_gene_name', values=external_gene_name[1]),
           error=function(x) data.frame(description='Could not query biomaRt!')) -> description 
  
  if(nrow(description)==0) # no description found
    return('No description of chosen feature found in biomaRt')
  
  description %>% unlist() %>% str_remove_all(pattern=' \\[Source.+$') %>% str_trunc(width=100)
}
