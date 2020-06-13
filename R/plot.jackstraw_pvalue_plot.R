#' Make an JackStraw pvalue plot of a reduction
#' 
#' Creates a plot of points and thresholds
#' 
#' @param id unique name of the element
#' 
#' @rdname jackstraw_pvalue_plot
#' 
jackstraw_pvalue_plot.ui <- function(id) {
  sprintf(fmt='### %s-jackstraw_pvalue_plot.ui', id) %>% message()

  module <- 'jackstraw_pvalue_plot'

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
  plotOutput(outputId=ns(id='jackstraw_pvalue_plot')) %>% withSpinner())
}

#' Produce the ggplot object for a JackStraw pvalue plot
#' 
#' @importFrom Seurat JS
#' 
#' @rdname jackstraw_pvalue_plot
#'
jackstraw_pvalue_plot.server <- function(input, output, session) {
  session$ns('') %>% sprintf(fmt='### %sjackstraw_pvalue_plot.server') %>% message()

  # get environments containing variables to run/configure this object
  collect_environments(id=parent.frame()$id, module='jackstraw_pvalue_plot') # provides `seuratvis_env`, `server_env` and `module_env`
  id <- parent.frame()$id
  ns <- NS(namespace='jackstraw_pvalue_plot')

  # render the elbow plot
  renderPlot(expr={
    # make sure these elements are defined
    req(seurat_object.reactions$seurat)
    req(input$reduction_method_picker)

    # send a message
    session$ns('') %>% sprintf(fmt='### %sjackstraw_pvalue_plot.server-renderPlot') %>% message('')
   
    # make variables for shorthand    
    seurat <- seurat_object.reactions$seurat
    reduction_name <- input$reduction_method_picker

    # make a plot
    as.data.frame(Seurat::JS(object=seurat[[reduction_name]], slot='overall')) %>% 
      ggplot() +
      aes(x=PC, y=-log10(Score)) +
      labs(x='Principle component', y='-log10(score)') +
      geom_smooth() +
      geom_point(shape=4) +
      theme_bw() +
      theme(legend.background=element_blank(),
            legend.justification=c(1,1),
            legend.position=c(1,1))}) -> output$jackstraw_pvalue_plot
}
