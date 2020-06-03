#' Make a boxplot of a feature
#' 
#' Creates a plot of points and thresholds
#' 
#' @param id unique name of the element
#' @param feature name of the column in the Seurat object to plot
#' 
#' @examples
#' 
#' \dontrun{
#' ui <- fluidPage(boxplot_plot.ui(id='page_name', feature='column_name'))
#' server <- function(input, output, session) {
#'   callModule(boxplot_plot.server, id='page_name-column_name')
#' }} 
#' 
#' @rdname boxplot_plot
#' 
boxplot_plot.ui <- function(id, feature) {
  sprintf(fmt='### %s-boxplot_plot.ui [%s]', id, feature) %>% message()

  id %<>% NS(id=feature) # combine the id and feature to allow multiple knee plots per id
  module <- 'boxplot_plot'

  # make unique id for this object
  ns <- NS(namespace=id)
  module_ns <- ns(id=module)

  # create an environment in the seuratvis namespace
  e <- new.env()
  e$id <- id
  e$feature <- feature # so we can plot the correct plot in the server
  assign(x=module_ns, val=e, envir=module_environments)

  # record the server(s) to call
  get0(env=module_servers_to_call, x=id) %>% append(sprintf(fmt='%s.server', module)) %>% assign(env=module_servers_to_call, x=id)

  # return ui element(s)
  plotOutput(outputId=ns(id='boxplot_plot')) %>% withSpinner()
}

#' Produce the ggplot object for a boxplot
#' 
#' @rdname boxplot_plot
#'
boxplot_plot.server <- function(input, output, session) {
  session$ns('') %>% sprintf(fmt='### %sboxplot_plot.server') %>% message()

  # get environments containing variables to run/configure this object
  collect_environments(id=parent.frame()$id, module='boxplot_plot') # provides `seuratvis_env`, `server_env` and `module_env`
  id <- parent.frame()$id

  # render the boxplot
  renderPlot(expr={
    sprintf(fmt='### boxplot_plot.server-renderPlot [%s]', id) %>% message()
   
    # get feature-specific plotting elements
    feature_plot <- NULL
    if(module_env$feature=='nCount_RNA') {
      # get thresholds of total UMI
      min_y <- filtering_parameters.reactions$total_umi_per_cell_min
      max_y <- filtering_parameters.reactions$total_umi_per_cell_max
      
      # start the boxplot
      seurat_object.reactions$n_umi_values %>%
        set_names('y') %>%
        ggplot() +
        aes(y=y) +
        labs(y='Total UMI per cell') +
        annotate(geom='rect', xmin=-Inf, xmax=Inf, ymin=min_y, ymax=max_y, fill='orange', alpha=0.2) +
        scale_y_continuous(labels=function(y) scales::comma(y, accuracy=1)) -> feature_plot
    } else if(module_env$feature=='nFeature_RNA') {
      # get thresholds of detected features
      min_y <- filtering_parameters.reactions$features_per_cell_min
      max_y <- filtering_parameters.reactions$features_per_cell_max
      
      # start the boxplot
      seurat_object.reactions$n_features_values %>%
        set_names('y') %>%
        ggplot() +
        aes(y=y) +
        labs(y='Features detected per cell') +
        annotate(geom='rect', xmin=-Inf, xmax=Inf, ymin=min_y, ymax=max_y, fill='orange', alpha=0.2) +
        scale_y_continuous(labels=function(y) scales::comma(y, accuracy=1)) -> feature_plot
    } else if(module_env$feature=='percent_mt') {
      # get max fraction Mt UMI
      max_y <- filtering_parameters.reactions$max_percent_mitochondria

      # start the boxplot
      seurat_object.reactions$proportion_mt_values %>%
        set_names('y') %>%
        ggplot() +
        aes(y=y) +
        labs(y='Proportion mitochondrial expression') +
        annotate(geom='rect', xmin=-Inf, xmax=Inf, ymin=0, ymax=max_y, fill='orange', alpha=0.2) -> feature_plot
    }

    # add the shared plotting elements
    feature_plot +
      aes(x=0) +
      geom_point(shape=16, size=0.6, colour='lightgrey', alpha=0.6, position=position_jitter(width=0.5)) +
      geom_boxplot(fill=NA, width=0.4, size=0.9, outlier.size=0.6) +
      coord_cartesian(xlim=c(-1, 1)) +
      theme_bw() +
      theme(axis.ticks.length.x=unit(0, 'pt'),
            axis.text.x=element_text(colour='white'),
            axis.title.x=element_text(colour='white'),
            panel.grid.major.x=element_blank(),
            panel.grid.minor.x=element_blank())}) -> output$boxplot_plot
}
