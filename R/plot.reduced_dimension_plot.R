#' Make a reduced dimension plot coloured by a feature
#' 
#' Creates a scatterplot of a reduced dimensionality map
#' 
#' @param id unique name of the element
#' @param feature either \code{selected_feature_resolution} or \code{picked_feature_values}
#' 
#' @examples
#' 
#' \dontrun{
#' ui <- fluidPage(reduced_dimension_plot.ui(id='page_name', feature='feature'))
#' server <- function(input, output, session) {
#'   callModule(reduced_dimension_plot.server, id='page_name-feature')
#' }} 
#' 
#' @rdname reduced_dimension_plo
#' 
reduced_dimension_plot.ui <- function(id, feature) {
  sprintf(fmt='### reduced_dimension_plot.ui [%s-%s]', id, feature) %>% message()

  id %<>% NS(id=feature) # combine the id and feature to allow multiple knee plots per id
  module <- 'reduced_dimension_plot'

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
  plotOutput(outputId=ns(id='reduced_dimension_plot')) %>% withSpinner()
}

#' Produce the ggplot object for a scatterplot
#' 
#' @rdname reduced_dimension_plo
#'
reduced_dimension_plot.server <- function(input, output, session) {
    sprintf('### reduced_dimension_plot.server') %>% message()

  # get environments containing variables to run/configure this object
  collect_environments(id=parent.frame()$id, module='reduced_dimension_plot') # provides `seuratvis_env`, `server_env` and `module_env`
  id <- parent.frame()$id
  session_server <- get(x='session', env=server_env)
  input_server <- get(x='input', env=server_env)

  # render the reduced dimension plot
  renderPlot(expr={
    sprintf(fmt='### reduced_dimension_plot.server-renderPlot [%s]', id) %>% message()
   
    # get feature-specific plotting elements
    output_plot <- NULL
    if(module_env$feature=='selected_cluster_resolution') {
      # get reduced dimension coordinates and currently selected cluster resolution and make the scatterplot
      cbind(seurat_object.reactions$dimred, seurat_object.reactions$picked_cluster_resolution_idents) %>%
        ggplot() +
        aes(x=DIMRED_1, y=DIMRED_2, colour=ident) +
        geom_hline(yintercept=0, colour='grey90') + geom_vline(xintercept=0, colour='grey90') +
        geom_point(size=seurat_object.reactions$point_size, alpha=seurat_object.reactions$opacity) +
        theme_void() +
        theme(legend.position='none') -> output_plot

      # if labels should be added, add them
      if(seurat_object.reactions$label_clusters) {
        cbind(seurat_object.reactions$dimred, seurat_object.reactions$picked_cluster_resolution_idents) %>%
          group_by(ident) %>%
          summarise(DIMRED_1=mean(DIMRED_1), DIMRED_2=mean(DIMRED_2)) -> data_labels

        output_plot +
          ggrepel::geom_label_repel(data=data_labels,
                                    mapping=aes(label=ident),
                                    colour='black',
                                    size=12/(14/5)) -> output_plot
      }
    } else if(module_env$feature=='picked_feature_values') {
      # get reduced dimension coordinates and the picked feature
      cbind(seurat_object.reactions$dimred, seurat_object.reactions$picked_feature_values) %>%
        rename(picked_feature_value=value) -> data

      # if the picked feature has numeric values
      if(is.numeric(data$picked_feature_value)) {
        # make the scatterplot
        data %>%
          arrange(picked_feature_value) %>%
          ggplot() +
          aes(x=DIMRED_1, y=DIMRED_2, colour=picked_feature_value) +
          geom_hline(yintercept=0, colour='grey90') + geom_vline(xintercept=0, colour='grey90') +
          geom_point(size=seurat_object.reactions$point_size, alpha=seurat_object.reactions$opacity) +
          scale_colour_gradient(low=input_server$`gene_highlighting-colour_palette-low`, high=input_server$`gene_highlighting-colour_palette-high`, limits=seurat_object.reactions$value_range_limits, oob=scales::squish) +
          theme_void() +
          theme(legend.position='none') -> output_plot
      } else { # the picked feature is factor-like
        # make the scatterplot and facet by picked feature value
        data %>%
          arrange(picked_feature_value) %>%
          ggplot() +
          aes(x=DIMRED_1, y=DIMRED_2, colour=picked_feature_value) +
          geom_point(size=seurat_object.reactions$point_size, alpha=seurat_object.reactions$opacity) +
          facet_wrap(~picked_feature_value, scales='free') +
          guides(colour=guide_legend(override.aes=list(size=3, shape=15))) +
          theme_void() +
          theme(legend.position='bottom', legend.title=element_blank()) -> output_plot
      }
    } else {
      sprintf('!!! cannot deal with %s', module_env$feature) %>% stop()
    }

    output_plot}) -> output$reduced_dimension_plot
}

# plotlyOutput('genes_highlighting-expression_per_cluster') %>% withSpinner()
# renderPlotly({plot_ly(z=~volcano) %>% add_surface()}) -> output$`genes_highlighting-expression_per_cluster`
