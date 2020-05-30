#' Display gene information
#' 
#' Displays gene of interest and biomaRt description
#' 
#' @param id unique name of the element
#' @param width integer value for width (rows sum to 12)
#' @param input,output,session used internally
#' 
#' @examples
#' 
#' \dontrun{
#' ui <- fluidPage(picked_feature_and_description_text_box.ui(id='page_name'))
#' server <- function(input, output, session) callModule(picked_feature_and_description_text_box.server, id='page_name')}
#' 
#' @rdname picked_feature_and_description_text_box
#' 
picked_feature_and_description_text_box.ui <- function(id, width=12) {
  module <- 'picked_feature_and_description_text_box'

  # record the server(s) to call
  get0(env=module_servers_to_call, x=id) %>% append(sprintf(fmt='%s.server', module)) %>% assign(env=module_servers_to_call, x=id)

  # return ui element(s)
  valueBoxOutput(outputId=NS(namespace=id, id='picked_feature_and_description'), width=width)
}

#' @rdname picked_feature_and_description_text_box
#' 
picked_feature_and_description_text_box.server <- function(input, output, session) {
  renderValueBox(env=parent.frame(n=2), quoted=FALSE, expr={
    message('### picked_feature_and_description_text_box.server')

    # create variables for shorthand
    picked_feature <- seurat_object.reactions$picked_feature

    # get gene description from biomaRt
    Misc(seurat_object.reactions$seurat, slot='mart') %>%
      get_formatted_gene_description(external_gene_name=picked_feature) -> subtitle

    # create output object
    list(value=picked_feature,
         subtitle=subtitle,
         icon=icon('jedi-order')) %>%
      modifyList(x=text_box_defaults()) %>%
      do.call(what=valueBox)}) -> output$picked_feature_and_description
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
  message('### get_formatted_gene_description')
  tryCatch(expr=getBM(mart=mart, attributes='description', filter='external_gene_name', values=external_gene_name[1]),
           error=function(x) data.frame(description='Could not query biomaRt!')) -> description 
  
  if(nrow(description)==0) # no description found
    return('No description of chosen feature found in biomaRt')
  
  description %>% unlist() %>% str_remove_all(pattern=' \\[Source.+$') %>% str_trunc(width=100)
}
