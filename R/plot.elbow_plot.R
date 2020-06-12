#' Make an elbow plot of a reduction
#' 
#' Creates a plot of points and thresholds
#' 
#' @param id unique name of the element
#' 
#' @rdname elbow_plot
#' 
elbow_plot.ui <- function(id) {
  sprintf(fmt='### %s-elbow_plot.ui', id) %>% message()

  module <- 'elbow_plot'

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
  tagList(plotOutput(outputId=ns(id='elbow_plot')) %>% withSpinner())
}

#' Produce the ggplot object for a elbow plot
#' 
#' @rdname elbow_plot
#'
elbow_plot.server <- function(input, output, session) {
  session$ns('') %>% sprintf(fmt='### %selbow_plot.server') %>% message()

  # get environments containing variables to run/configure this object
  collect_environments(id=parent.frame()$id, module='elbow_plot') # provides `seuratvis_env`, `server_env` and `module_env`
  id <- parent.frame()$id
  ns <- NS(namespace='elbow_plot')

  # render the elbow plot
  renderPlot(expr={
    # make sure these elements are defined
    req(seurat_object.reactions$seurat)

    # send a message
    session$ns('') %>% sprintf(fmt='### %selbow_plot.server-renderPlot') %>% message('')
   
    # make variables for shorthand    
    seurat <- seurat_object.reactions$seurat
    reduction_name <- input$reduction_method_picker

    data.frame(project_name=Project(seurat),
               Y=Stdev(object=seurat, reduction=reduction_name)) %>%
      mutate(X=seq(n())) -> data

    stdev <- Stdev(object=seurat, reduction=reduction_name)
    pct <- stdev / sum(stdev) * 100 # Determine percent of variation associated with each PC
    cumu <- cumsum(pct) # Calculate cumulative percents for each PC
    
    co1 <- which(cumu > 90 & pct < 5)[1] # Determine which PC exhibits cumulative percent greater than 90% and % variation associated with the PC as less than 5
    co2 <- sort(which((pct[1:length(pct) - 1] - pct[2:length(pct)]) > 0.1), decreasing=TRUE)[1] + 1 # Determine the difference between variation of PC and subsequent PC

    data.frame(co1=co1, co2=co2) %>%
      mutate(pcs=pmin(co1, co2),
             min=pmin(co1, co2),
             max=pmax(co1, co2)) -> dimensions

    data %>%
      mutate(selected=cut(x=X, breaks=c(0, co1, co2, length(stdev)+1))) %>%
      ggplot(data=.) +
      aes(x=X, y=Y, colour=selected) +
      labs(x='Principal component', y='Standard deviation', colour='Selected') +
      geom_line(colour='grey', size=1) +
      geom_point(size=2) +
      theme_bw() +
      theme(aspect.ratio=1,
            legend.position='bottom',
            legend.title=element_blank(),
            panel.grid.minor=element_blank(),
            strip.background=element_rect(fill=NA))}) -> output$elbow_plot
}
