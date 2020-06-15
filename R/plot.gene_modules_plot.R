#' Make a gene module score density plot
#' 
#' Creates a density plot of a gene module score
#' 
#' @param id unique name of the element
#' 
#' @rdname gene_module_score_plot
#' 
gene_module_score_plot.ui <- function(id) {
  sprintf(fmt='### %s-elbow_plot.ui', id) %>% message()

  module <- 'gene_module_score_plot'

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
  tagList(plotOutput(outputId=ns(id='gene_module_score_ridges')) %>% withSpinner())
}

#' Produce the ggplot object for a gene module score
#' 
#' @import ggridges
#' 
#' @rdname gene_module_score_plot
#'
gene_module_score_plot.server <- function(input, output, session, seurat, ...) {
  session$ns('') %>% sprintf(fmt='### %selbow_plot.server') %>% message()

  # get environments containing variables to run/configure this object
  collect_environments(id=parent.frame()$id, module='gene_module_score_plot') # provides `seuratvis_env`, `server_env` and `module_env`
  tab <- parent.frame()$id
  ns <- NS(namespace='gene_module_score_plot')

  updateRadioGroupButtons(session=session, inputId='feature_type',
                          choices=list(`Gene modules`='gene_modules'),
                          selected='gene_modules')

  # render the elbow plot
  renderPlot(expr={
    # make sure these elements are defined
    req(seurat$object)

    # send a message
    session$ns('') %>% sprintf(fmt='### %sgene_module_score_plot.server-renderPlot') %>% message()
   
    # make variables for shorthand    
    object <- seurat$object
    reduction_name <- input$reduction_method_picker

    # plot module (y) and score (x)
    # cbind(seurat$picked_feature_values[[tab]] %>% set_names('score'),
    #       seurat$picked_cluster_resolution_idents) %>%
    #   gather(key=gene_module, value=module_score, -cluster_id) %>%
    #   ggplot() +
    #   aes(x=module_score, y=gene_module, fill=gene_module) +
    #   labs(title=str_c(Project(seurat), cluster_set, cluster_id, sep=' / ')) +
    #   geom_density_ridges(alpha=0.8, size=0.7, colour='black') +
    #   scale_fill_cyclical(values=c('grey50','grey70')) +
    #   theme_bw() +
    #   theme(axis.text.y=element_text(vjust=0), panel.grid.major.x=element_blank(), panel.grid.minor.x=element_blank())

    # plot cluster (y) and score(x)
    cbind(seurat$picked_cluster_resolution_idents,
          seurat$picked_feature_values[[tab]]) %>%
      set_names(c('cluster_id', 'module_score')) %>%
      ggplot() +
      aes(x=module_score, y=cluster_id, fill=cluster_id) +
      labs(title=str_c(seurat$formatted_project, input$cluster_resolution_picker, input$picked_feature, sep=' / ')) +
      geom_density_ridges(alpha=0.8, size=0.7, colour='black') +
      scale_fill_cyclical(values=c('grey50','grey70')) +
      theme_bw() +
      theme(axis.text.y=element_text(vjust=0), panel.grid.major.x=element_blank(), panel.grid.minor.x=element_blank())



######
# need to set cluster idents to one
# need to fetch multiple features
# plot the reults
######




    }, height=800) -> output$gene_module_score_ridges
}
