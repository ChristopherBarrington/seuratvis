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
  sprintf(fmt='### %s-reduced_dimension_plot.ui [%s]', id, feature) %>% message()

  id %<>% NS(id=feature) # combine the id and feature to allow multiple UMAP plots per id
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
#' @import ggrepel
#' 
#' @rdname reduced_dimension_plot
#'
reduced_dimension_plot.server <- function(input, output, session) {
  session$ns('') %>% sprintf(fmt='### %sreduced_dimension_plot.server') %>% message()

  # get environments containing variables to run/configure this object
  collect_environments(id=parent.frame()$id, module='reduced_dimension_plot') # provides `seuratvis_env`, `server_env` and `module_env`
  id <- parent.frame()$id
  session_server <- get(x='session', env=server_env)
  input_server <- get(x='input', env=server_env)

  # render the reduced dimension plot
  renderPlot(expr={

    req(selections.rv[[session$ns('dimred') %>% parse_ns_label()]])

    # send a message
    session$ns('') %>% sprintf(fmt='### %sreduced_dimension_plot.server-renderPlot') %>% message()

    # collect args from selections.rv
    c('dimred', 'picked_feature_values') %>%
      purrr::set_names() %>%
      sapply(session$ns) %>%
      sapply(parse_ns_label) %>%
      lapply(function(x) selections.rv[[x]]) -> args

    # because multiple plots are using the same `input` tag...
    #! TODO: could make each umap type a different module?? ie differnet ui id tags? in different servers perhaps?
    args$opacity <- input_server[[session$ns('opacity_slider') %>% parse_ns_label()]]
    args$point_size <- input_server[[session$ns('point_size_slider') %>% parse_ns_label()]]
    args$value_range_limits <- input_server[[session$ns('value_range') %>% parse_ns_label()]]
    args$cluster_id_picker <- input_server[[session$ns('cluster_id_picker') %>% parse_ns_label()]]
    args$label_clusters <- input_server[[session$ns('label_clusters') %>% parse_ns_label()]]

    # make a base plot
    cbind(args$dimred,
          seurat_object.reactions$picked_cluster_resolution_idents,
          {args$picked_feature_values %>% rename(picked_feature_value=value)}) %>%
      mutate(is_selected_cluster_id=ident %in% args$cluster_id_picker) %>%
      arrange(picked_feature_value) %>%
      ggplot() +
      aes(x=DIMRED_1, y=DIMRED_2) +
      geom_hline(yintercept=0, colour='grey90') + geom_vline(xintercept=0, colour='grey90') +
      geom_point(size=args$point_size, alpha=args$opacity) +
      theme_void() +
      theme(legend.position='none', legend.text=element_blank()) -> output_plot

    # get feature-specific plotting elements
    include_legend <- FALSE
    if(module_env$feature=='selected_cluster_resolution') {
      # get reduced dimension coordinates and currently selected cluster resolution and make the scatterplot
      output_plot +
        aes(colour=ident) -> output_plot

      # if labels should be added, add them
      if(args$label_clusters) {
        output_plot$data %>%
          group_by(ident) %>%
          summarise(DIMRED_1=mean(DIMRED_1), DIMRED_2=mean(DIMRED_2)) -> data_labels

        output_plot +
          geom_label_repel(data=data_labels,
                           mapping=aes(label=ident),
                           colour='black',
                           size=12/(14/5)) -> output_plot
      }
      include_legend <- FALSE
    } else if(module_env$feature=='picked_feature_values') {
      # if the picked feature has numeric values
      if(is.numeric(output_plot$data$picked_feature_value)) {
        # make the scatterplot
        output_plot +
          aes(colour=picked_feature_value) -> output_plot

        c_min <- session$ns('low') %>% parse_ns_label() %>% pluck(.x=plotting_options.rv$colours) # TODO: this is dependent on the label names!
        c_mid <- 'white'
        c_max <- session$ns('high') %>% parse_ns_label() %>% pluck(.x=plotting_options.rv$colours) # TODO: this is dependent on the label names!
        c_range_limits <- args$value_range_limits

        colour_gradient <- scale_colour_gradient(low=c_min, high=c_max, limits=c_range_limits, oob=scales::squish)
        if(c_range_limits %>% sign() %>% Reduce(f='*') %>% equals(-1)) {
          colour_gradient <- scale_colour_gradientn(colours=c(low=c_min, mid=c_mid, high=c_max), 
                                                    values={c_range_limits %>% c(0) %>% sort() %>% scales::rescale()},
                                                    limits=c_range_limits, breaks=0)
          output_plot <- output_plot + theme(legend.text=element_text())
        }

        output_plot <- output_plot + colour_gradient
        include_legend <- TRUE
      } else { # the picked feature is factor-like
        # make the scatterplot and facet by picked feature value
        output_plot +
          aes(colour=picked_feature_value) +
          facet_wrap(~picked_feature_value, scales='fixed') +
          theme(panel.background=element_rect(fill=NA, colour='black'),
                strip.background=element_rect(fill=NA, colour='black'),
                legend.position='none') -> output_plot
        include_legend <- FALSE
      }
    } else if(module_env$feature=='selected_cluster_ids') {
      output_plot +
        aes(colour=is_selected_cluster_id, alpha=is_selected_cluster_id) +
        scale_colour_manual(values=c(`FALSE`='grey90', `TRUE`='darkorange1')) +
        scale_alpha_manual(values=c(`FALSE`=0.1, `TRUE`=1)) -> output_plot
      output_plot$data %<>% arrange(is_selected_cluster_id)
      output_plot$layers[[3]]$aes_params$alpha <- NULL
      include_legend <- FALSE
    } else {
      sprintf('!!! cannot deal with %s', module_env$feature) %>% stop()
    }

    if(include_legend)
      output_plot +
         theme(legend.justification=c(1,0),
               legend.position=c(1,0),
               legend.direction='horizontal',
               legend.title=element_blank()) -> output_plot

    output_plot}) -> output$reduced_dimension_plot
}

# plotlyOutput('genes_highlighting-expression_per_cluster') %>% withSpinner()
# renderPlotly({plot_ly(z=~volcano) %>% add_surface()}) -> output$`genes_highlighting-expression_per_cluster`
