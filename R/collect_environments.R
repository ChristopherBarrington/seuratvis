#' Collect environments for modules
#' 
#' Creates environments for module-specific and seuratvis-wide variables
#' 
#' @param module name of the module used to determine module-specific environment to load
#' 
#' @return none
#' 
#' @details assigns a \code{seuratvis_env} and \code{module_env} environment in the calling function, including copying the \code{id} variable that is passed to \code{callModule()}
#' 
collect_environments <- function(module) {
  # get the environment using the shiny object identifiers (namespace and id)
  callmodule_env <- parent.frame(n=2)
  callmodule_env %>%
    pluck('id') %>%
    NS(id=module) %>%
    sprintf(fmt='%s.env') %>%
    getFromNamespace(ns='seuratvis') -> module_env

  # assign environments back to the calling server function
  assign(x='seuratvis_env', value=parent.frame(n=3), envir=parent.frame(n=1))
  assign(x='module_env', value=module_env, envir=parent.frame(n=1))

  invisible('done')
}
