#' Collect environments for modules
#' 
#' Creates environments for module-specific and seuratvis-wide variables
#' 
#' @param module name of the module used to determine module-specific environment to load
#' @param id name of the UI element
#' 
#' @return none
#' 
#' @details assigns a \code{seuratvis_env} and \code{module_env} environment in the calling function, including copying the \code{id} variable that is passed to \code{callModule()}
#' 
collect_environments <- function(module, id) {
  NS(namespace=id, id=module) %>% get(envir=module_environments) -> module_env

  assign(x='seuratvis_env', value=parent.env(environment()), envir=parent.frame(n=1))
  assign(x='module_env', value=module_env, envir=parent.frame(n=1))
  assign(x='server_env', value=parent.frame(n=3), envir=parent.frame(n=1))

  invisible('done')
}
