#' Find Seurat objects
#' 
#' Look for Seurat objects in environments
#' 
#' @details The function uses \code{ls()} to get to first get a list of environments in the global environment. For each environment, \code{ls()} is again used to list the variables. If a variable is a \code{Seurat} object (its \code{class} is \code{Seurat}), it is retained in the output \code{data.frame}.
#' 
#' @return A \code{data.frame}
#' 
find_seurat_objects <- function() {
  ls(envir=globalenv()) %>%
    sapply(function(O) get(x=O, envir=globalenv()) %>% class()) %>%
    enframe() %>%
    plyr::dlply(~value, pluck, 'name') %>%
    pluck('environment') %>%
    rev() %>%
    sapply(get, envir=globalenv()) %>%
    append(list(`globalenv()`=globalenv())) %>%
    rev() %>%
    lapply(function(E) {ls(envir=E) %>% sapply(function(O) get(x=O, envir=E) %>% class()) %>% enframe() %>% plyr::dlply(~value, pluck, 'name') %>% pluck('Seurat')}) %>%
    plyr::ldply(.id='env', enframe) %>%
    dplyr::select(-name) %>%
    unite(col='choiceValue', sep='$', env, value, remove=FALSE) %>%
    (function(x) {
      if(length(unique(x$env))==1) {
        x %>% mutate(choiceName=value)
      } else {
        x %>% mutate(choiceName=str_c(str_remove_all(string=env, pattern='\\(\\)$'), value, sep=' : '))
      }}) %>%
    arrange(choiceName)
}
