
#' Get description of a gene
#' 
#' Use biomaRt to get formatted description of a gene from gene name
#' 
#' @param external_gene_name name of gene to query
#' @param mart a \code{biomaRt} object
#' 
#' @return a character string of the \code{attribute} from the \code{mart} for the \code{external_gene_name}.
#' 
get_formatted_gene_description <- function(external_gene_name, mart) {
  description <- getBM(mart=mart, attributes='description', filter='external_gene_name', values=external_gene_name[1])
  
  if(nrow(description)==0) # no description found
    return('Selected gene')
  
  str_remove_all(string=description, pattern=' \\[Source.+$') %>% deframe() %>% str_trunc(width=100)
}
