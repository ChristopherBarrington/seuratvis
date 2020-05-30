#' Make a density plot of a feature
#' 
#' Creates a density plot and thresholds
#' 
#' @param id unique name of the element
#' @param feature name of the column in the Seurat object to plot
#' 
#' @examples
#' 
#' \dontrun{
#' ui <- fluidPage(density_plot.ui(id='page_name', feature='column_name'))
#' server <- function(input, output, session) {
#'   callModule(density_plot.server, id='page_name-column_name')
#' }} 
#' 
#' @rdname density_plot
#' 
density_plot.ui <- function(id, feature) {
  sprintf(fmt='### density_plot.ui [%s-%s]', id, feature) %>% message()

  id %<>% NS(id=feature) # combine the id and feature to allow multiple knee plots per id
  module <- 'density_plot'

  # make unique id for this object
  ns <- NS(namespace=id)
  module_ns <- ns(id=module)

  # create an environment in the seuratvis namespace
  e <- new.env()
  e$id <- id
  e$feature <- feature
  assign(x=module_ns, val=e, envir=module_environments)

  # record the server(s) to call
  get0(env=module_servers_to_call, x=id) %>% append(sprintf(fmt='%s.server', module)) %>% assign(env=module_servers_to_call, x=id)

  # return ui element(s)
  plotOutput(outputId=ns(id='density_plot'), brush=brushOpts(id=ns('brush'), direction='x')) %>% withSpinner()
}

#' Produce the ggplot object for a boxplot
#' 
#' @rdname density_plot
#'
density_plot.server <- function(input, output, session) {
  sprintf('### density_plot.server') %>% message()

  # get environments containing variables to run/configure this object
  collect_environments(id=parent.frame()$id, module='density_plot') # provides `seuratvis_env`, `server_env` and `module_env`
  id <- parent.frame()$id
  session_server <- get(x='session', env=server_env)

  # render the density plot
  renderPlot(expr={
    sprintf(fmt='### density_plot.server-renderPlot [%s]', id) %>% message()
   
    # get feature-specific plotting elements
    feature_plot <- NULL
    if(module_env$feature=='nCount_RNA') {
      # start the density plot
      FetchData(seurat_object.reactions$seurat, 'nCount_RNA') %>%
        set_names('y') %>%
        ggplot() +
        aes(x=y) +
        labs(x='Total UMI per cell', y='Density') +
        stat_density(geom='line', trim=TRUE, size=2) +
        scale_x_log10(breaks=major_breaks_log10, minor_breaks=minor_breaks_log10, labels=function(x) scales::comma(x, accuracy=1))+
        theme_bw()+
        theme() -> feature_plot
    } else if(module_env$feature=='nFeature_RNA') {
      # start the density plot
      FetchData(seurat_object.reactions$seurat, 'nFeature_RNA') %>%
        set_names('y') %>%
        ggplot() +
        aes(x=y) +
        labs(x='Features detected per cell', y='Density') +
        stat_density(geom='line', trim=TRUE, size=2) +
        scale_x_log10(breaks=major_breaks_log10, minor_breaks=minor_breaks_log10, labels=function(x) scales::comma(x, accuracy=1))+
        theme_bw()+
        theme() -> feature_plot
    } else if(module_env$feature=='percent_mt') {
      # start the density plot
      seurat_object.reactions$percent_mt %>%
        set_names('y') %>%
        ggplot() +
        aes(x=y) +
        labs(x='Proportion mitochondrial expression', y='Density') +
        stat_density(geom='line', trim=TRUE, size=2) +
        theme_bw()+
        theme() -> feature_plot
    }

    # add the shared plotting elements
    feature_plot}) -> output$density_plot

  # react to the brush
  observeEvent(eventExpr=input$brush, handlerExpr={
    sprintf(fmt='### density_plot.server-observeEvent-input$brush [%s]', id) %>% message()

    # update feature-specific reactives and ui elements
    if(module_env$feature=='nFeature_RNA') {
      # get brush positions, round up or down
      low <- floor(input$brush$xmin)
      high <- ceiling(input$brush$xmax)

      # update the ui
      for(id in module_environments$features_per_cell_filters$id) {
        updateTextInput(session=session_server, inputId=NS(id, 'min_features'), value=low)
        updateTextInput(session=session_server, inputId=NS(id, 'max_features'), value=high)
      }

      # update the reactive
      filtering_parameters.reactions$features_per_cell_min <- low
      filtering_parameters.reactions$features_per_cell_max <- high
    } else if(module_env$feature=='nCount_RNA') {
      # get brush positions, round up or down
      low <- floor(input$brush$xmin)
      high <- ceiling(input$brush$xmax)

      # update the ui
      for(id in module_environments$total_umi_per_cell_filters$id) {
        updateTextInput(session=session_server, inputId=NS(namespace=id, id='min_umis'), value=low)
        updateTextInput(session=session_server, inputId=NS(namespace=id, id='max_umis'), value=high)
      }

      # update the reactive
      filtering_parameters.reactions$total_umi_per_cell_min <- low
      filtering_parameters.reactions$total_umi_per_cell_max <- high
    } else if(module_env$feature=='percent_mt') {
      # get the max position of the brush to 1dp
      high <- round(input$brush$xmax+0.05, digits=1)

      # update the ui
      for(id in module_environments$percent_mt_per_cell_filters$id)
        updateTextInput(session=session_server, inputId=NS(namespace=id, id='max_percent_mt'), value=high)

      # update the reactive
      filtering_parameters.reactions$percent_mt_per_cell_max <- high
    }})
}
