#' Make a knee plot of a feature
#' 
#' Creates a plot of points and thresholds
#' 
#' @param id unique name of the element
#' @param feature name of the column in the Seurat object to plot
#' 
#' @examples
#' 
#' \dontrun{
#' ui <- fluidPage(knee_plot.ui(id='page_name', feature='column_name'))
#' server <- function(input, output, session) {
#'   callModule(knee_plot.server, id='page_name-column_name')
#' }} 
#' 
#' @rdname knee_plot
#' 
knee_plot.ui <- function(id, feature) {
  sprintf(fmt='### knee_plot.ui [%s-%s]', id, feature) %>% message()

  id %<>% NS(id=feature) # combine the id and feature to allow multiple knee plots per id
  module <- 'knee_plot'

  # make unique id for this object
  ns <- NS(namespace=id)
  module_ns <- ns(id=module)

  # create an environment in the seuratvis namespace
  e <- new.env()
  e$id <- id
  e$feature <- feature
  assign(x=module_ns, val=e, envir=module_environments)

  module_environments$knee_plots$ns %<>% c(module_ns)
  module_environments$knee_plots$id %<>% c(id)

  # return ui element(s)
  plotOutput(outputId=ns(id='knee_plot')) %>% withSpinner()
}

#' Produce the ggplot object for a knee plot
#' 
#' @rdname knee_plot
#'
knee_plot.server <- function(input, output, session) {
  message('### knee_plot.server')

  # get environments containing variables to run/configure this object
  collect_environments(id=parent.frame()$id, module='knee_plot') # provides `seuratvis_env`, `server_env` and `module_env`
  id <- parent.frame()$id

  # render the knee plot
  renderPlot(expr={
    sprintf(fmt='### knee_plot.server-renderPlot [%s]', id) %>% message()
   
    # get feature-specific plotting elements
    feature_plot <- NULL
    if(module_env$feature=='nCount_RNA') {
      # get thresholds of total UMI
      min_value <- filtering_parameters.reactions$total_umi_per_cell_min
      max_value <- filtering_parameters.reactions$total_umi_per_cell_max
      
      # start the knee plot
      FetchData(seurat_object.reactions$seurat, 'nCount_RNA') %>%
        set_names('y') %>%
        arrange(desc(y)) %>%
        mutate(x=seq(n()),
               pass=between(x=y, left=min_value, right=max_value)) %>%
        ggplot() +
        aes(x=x, y=y, colour=pass) +
        labs(y='Total UMI per cell') +
        geom_hline(yintercept=min_value, size=1) +
        scale_y_log10(breaks=major_breaks_log10, minor_breaks=minor_breaks_log10, labels=function(y) scales::comma(y, accuracy=1)) -> feature_plot
    } else if(module_env$feature=='nFeature_RNA') {
      # get thresholds of detected features
      min_value <- filtering_parameters.reactions$features_per_cell_min
      max_value <- filtering_parameters.reactions$features_per_cell_max
      
      # start the knee plot
      FetchData(seurat_object.reactions$seurat, 'nFeature_RNA') %>%
        set_names('y') %>%
        arrange(desc(y)) %>%
        mutate(x=seq(n()),
               pass=between(x=y, left=min_value, right=max_value)) %>%
        ggplot() +
        aes(x=x, y=y, colour=pass) +
        labs(y='Features detected per cell') +
        geom_hline(yintercept=min_value, size=1) +
        scale_y_log10(breaks=major_breaks_log10, minor_breaks=minor_breaks_log10, labels=function(y) scales::comma(y, accuracy=1)) -> feature_plot
    } else if(module_env$feature=='percent_mt') {
      # get max fraction Mt UMI
      max_value <- filtering_parameters.reactions$max_percent_mitochondria

      # start the knee plot
      seurat_object.reactions$percent_mt %>%
        set_names('y') %>%
        arrange(desc(y)) %>%
        mutate(x=seq(n()),
               pass=y<=max_value) %>%
        ggplot()+
        aes(x=x, y=y, colour=pass)+
        labs(y='Proportion mitochondrial expression') -> feature_plot
    }

    # add the shared plotting elements
    feature_plot +
      labs(x='Ranked cells', colour='Cells passing threshold') +
      geom_hline(yintercept=max_value, size=1) +
      geom_point(alpha=1, shape=16) +
      scale_x_log10(breaks=major_breaks_log10, minor_breaks=minor_breaks_log10, labels=function(x) scales::comma(x, accuracy=1)) +
      scale_colour_brewer(palette='Set1', direction=1) +
      theme_bw() +
      theme(legend.position='none', legend.title=element_blank())}) -> output$knee_plot
}
