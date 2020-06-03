#' Make a plot of a feature values in clusters
#' 
#' Creates a plot to show median and range of selected feature values in cluster idents
#' 
#' @param id unique name of the element
#' 
#' @examples
#' 
#' \dontrun{
#' ui <- fluidPage(feature_values_per_cluster.ui(id='page_name'))
#' server <- function(input, output, session) {
#'   callModule(feature_values_per_cluster.server, id='page_name')
#' }} 
#' 
#' @rdname feature_values_per_cluster_plot
#' 
feature_values_per_cluster_plot.ui <- function(id) {
  sprintf(fmt='### %s-feature_values_per_cluster_plot.ui', id) %>% message()

  module <- 'feature_values_per_cluster_plot'

  # make unique id for this object
  ns <- NS(namespace=id)
  module_ns <- ns(id=module)

  # create an environment in the seuratvis namespace
  e <- new.env()
  e$id <- id
  assign(x=module_ns, val=e, envir=module_environments)

  # record the server(s) to call
  get0(env=module_servers_to_call, x=id) %>% append(sprintf(fmt='%s.server', module)) %>% assign(env=module_servers_to_call, x=id)

  # return ui element(s)
  plotOutput(outputId=ns(id='feature_values_per_cluster_plot')) %>% withSpinner()
}

#' Produce the ggplot object for feature values in clusters plot
#' 
#' @rdname feature_values_per_cluster_plot
#'
feature_values_per_cluster_plot.server <- function(input, output, session) {
  session$ns('') %>% sprintf(fmt='### %sfeature_values_per_cluster_plot.server') %>% message()

  # get environments containing variables to run/configure this object
  collect_environments(id=parent.frame()$id, module='feature_values_per_cluster_plot') # provides `seuratvis_env`, `server_env` and `module_env`
  id <- parent.frame()$id

  # render the knee plot
  renderPlot(expr={
    # send a message
    session$ns('') %>% sprintf(fmt='### %sfeature_values_per_cluster_plot.server-renderPlot') %>% message('')
   
    # get the data to plot
    cbind(seurat_object.reactions$picked_cluster_resolution_idents,
          selections.rv[[session$ns('picked_feature_values')]]) %>%
      mutate(x={as.character(ident) %>% as.numeric()}) -> data

    # if the picked feature has numeric values
    if(is.numeric(data$value)) {
      # make a values range-type plot
      data %>%
        # filter(value>0) %>%
        group_by(ident, x) %>%
        summarise(q25=quantile(value, 0.25), q75=quantile(value, 0.75), median=median(value)) %>%
        mutate(iqr=q75-q25, lower=q25-1.5*iqr, upper=q75+1.5*iqr) -> cluster_data_summary

      cluster_data_summary %>%
        gather(key=key, value=y, lower, upper) %>%
        ggplot() +
        aes(x=x, y=y, colour=ident) +
        labs(y='Feature value (median ± 1.5x IQR)') +
        geom_line(size=1) +
        geom_point(mapping=aes(y=median), colour='black', shape=20, size=3) -> feature_plot
    } else { # the picked feature is factor-like
      # make a frequency bar plot
      data %>%
        ggplot() +
        aes(x=x, fill=value) +
        labs(y='Frequency') +
        geom_bar(position='dodge') -> feature_plot
    }

    # add the shared plotting elements
    feature_plot +
        labs(x='Cluster identifier') +
        theme_bw() +
        theme(legend.position='none', panel.grid.major.x=element_blank(), panel.grid.minor.x=element_blank())}) -> output$feature_values_per_cluster_plot
}
