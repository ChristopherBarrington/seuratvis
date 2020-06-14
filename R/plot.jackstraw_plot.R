#' Make a JackStraw plot of a reduction
#' 
#' Creates a JackStraw plot
#' 
#' @param id unique name of the element
#' 
#' @rdname jackstraw_plot
#' 
jackstraw_plot.ui <- function(id) {
  sprintf(fmt='### %s-jackstraw_plot.ui', id) %>% message()

  module <- 'jackstraw_plot'

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
  tagList(plotOutput(outputId=ns(id='jackstraw_plot')) %>% withSpinner())
}

#' Produce the ggplot object for a JackStraw plot
#' 
#' @importFrom Seurat JS
#' 
#' @rdname jackstraw_plot
#'
jackstraw_plot.server <- function(input, output, session, ...) {
  session$ns('') %>% sprintf(fmt='### %sjackstraw_plot.server') %>% message()

  # get environments containing variables to run/configure this object
  collect_environments(id=parent.frame()$id, module='jackstraw_plot') # provides `seuratvis_env`, `server_env` and `module_env`
  id <- parent.frame()$id
  ns <- NS(namespace='jackstraw_plot')

  # render the elbow plot
  renderPlot(expr={
    # make sure these elements are defined
    req(seurat_object.reactions$seurat)

    # send a message
    session$ns('') %>% sprintf(fmt='### %sjackstraw_plot.server-renderPlot') %>% message('')
   
    # make variables for shorthand    
    seurat <- seurat_object.reactions$seurat
    reduction_name <- input$reduction_method_picker
    components_range <- input$principle_components_slider
    min_component <- min(components_range)
    max_component <- max(components_range)

    Seurat::JS(object=seurat[[reduction_name]], slot='empirical') %>%
      as.data.frame() %>%
      rownames_to_column('Contig') %>%
      gather(key='PC', value='Value', -Contig) %>%
      mutate(PC={stringr::str_remove(PC, '^PC') %>% as.numeric()}) %>%
      filter(dplyr::between(PC, left=min_component, right=max_component)) %>%
      left_join(y=as.data.frame(Seurat::JS(object=seurat[[reduction_name]], slot='overall')), by='PC') %>%
      mutate(PC_colour=sprintf('PC %d: %1.3g', PC, Score),
             project_name=Project(seurat)) %>%
      mutate(PC_colour=factor(PC_colour, levels={unique(PC_colour) %>% str_sort(numeric=TRUE)})) %>%
      ggplot() +
      aes(sample=Value, colour=PC_colour) +
      labs(x='Theoretical [runif(1000)]', y='Empirical', colour='PC pvalue') +
      stat_qq(distribution=qunif, alpha=0.5) +
      geom_abline(intercept=0, slope=1, linetype='dashed') +
      coord_flip() +
      theme_bw() +
      theme(legend.position='none')}) -> output$jackstraw_plot
}
