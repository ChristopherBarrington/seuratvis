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
    unlist() %>%
    enframe() %>%
    plyr::dlply(~value, pluck, 'name') %>%
    pluck('environment') %>% # select the environments from the RGlobalEnv
    rev() %>%
    sapply(get, envir=globalenv()) %>%
    append(list(`globalenv()`=globalenv())) %>% # make a list of all environments and RGlobalEnv
    rev() %>%
    lapply(function(E) {ls(envir=E) %>% sapply(function(O) get(x=O, envir=E) %>% class()) %>% unlist() %>% enframe() %>% plyr::dlply(~value, pluck, 'name') %>% pluck('Seurat')}) %>% # get the class of objects in the environments and select out the Seurat objects
    plyr::ldply(.id='env', enframe) %>%
    dplyr::select(-name) %>%
    unite(col='choiceValue', sep='$', env, value, remove=FALSE) %>%
    (function(x) {
      if(length(unique(x$env))==1) {
        x %>% mutate(choiceName=value)
      } else {
        x %>% mutate(choiceName=str_c(str_remove_all(string=env, pattern='\\(\\)$'), value, sep=' : '))
      }}) %>%
    arrange(choiceName) %>%
    mutate(env=as.character(env)) -> available_objects # add varaibles for where to find the objects and what to call them

  available_objects
}
