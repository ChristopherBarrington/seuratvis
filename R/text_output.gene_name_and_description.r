#' 
#' 
picked_feature_and_description_text_box.ui <- function(id, width=12) {
  valueBoxOutput(outputId=NS(id, 'picked_feature_and_description'), width=width)
}

#'
#' 
picked_feature_and_description_text_box.server <- function(input, output, session, seurat, picked_feature) {
  renderValueBox(expr={
    req(picked_feature$name)
    # req(seurat$mart)
    
    # create variables for shorthand
    picked_feature <- picked_feature$name

    # get gene description from biomaRt
    get_formatted_gene_description(external_gene_name=picked_feature, mart=seurat$mart) -> subtitle

    # create output object
    list(value=picked_feature,
         subtitle=subtitle,
         icon=icon('jedi-order')) %>%
      modifyList(x=list(color='purple')) %>%
      do.call(what=valueBox)}) -> output$picked_feature_and_description
}

#' 
#' 
get_formatted_gene_description <- function(external_gene_name, mart) {
  tryCatch(expr=getBM(mart=mart, attributes='description', filter='external_gene_name', values=external_gene_name[1]),
           error=function(x) data.frame(description='Could not query biomaRt!')) -> description 
  
  if(nrow(description)==0) # no description found
    return('No description of chosen feature found in biomaRt')
  
  description %>% unlist() %>% str_remove_all(pattern=' \\[Source.+$') %>% str_trunc(width=100)
}
