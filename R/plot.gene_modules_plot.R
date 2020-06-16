#' Make a gene module score density plot
#' 
#' Creates a density plot of a gene module score
#' 
#' @param id unique name of the element
#' 
#' @rdname gene_module_score_plot
#' 
gene_module_score_plot.ui <- function(id) {
  sprintf(fmt='### %s-gene_module_score_plot.ui', id) %>% message()

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

  # make the ui elements
  #! TODO: can this be generalised?
  ## scores in a cluster tab panel
  tabpanel_name <- str_c(id, 'score_in_clusters', sep='_')
  tabpanel_ns <- NS(namespace=tabpanel_name)
  tabPanel(title='Score in clusters', value='score_in_clusters',
    column(width=8, 
           div(style="display: inline-block;vertical-align:top; width: 150px;", cluster_id_picker.ui(id=tabpanel_name, opts=list(multiple=TRUE))),
           div(style="display: inline-block;vertical-align:top; width: 150px;", feature_picker.ui(id=tabpanel_name, gene_modules_opts=list(multiple=FALSE), include_feature_type=FALSE, include_values_range=FALSE)),
           plotOutput(outputId=tabpanel_ns('score_in_clusters_ridges')) %>% withSpinner()),
    column(width=4, boxPlus(title='Options', width=12, closable=FALSE,
                            cluster_resolution_picker.ui(id=tabpanel_name)),
                            # add_gene_module.ui(id=tab)),
                    boxPlus(title='Map overview', width=12, closable=FALSE,
                            reduced_dimension_plot.ui(id=tabpanel_name, feature='selected_cluster_ids'),
                            reduction_method_picker.ui(id=tabpanel_name)),
                    boxPlus(title='Gene modules', width=12, closable=FALSE,
                            show_genes_in_module.ui(id=tabpanel_name)))) -> score_in_clusters
  get0(env=module_servers_to_call, x=tabpanel_name) %>% append(sprintf(fmt='%s.server', module)) %>% assign(env=module_servers_to_call, x=tabpanel_name)

  ## score in clusters tab panel
  tabpanel_name <- str_c(id, 'scores_in_cluster', sep='_')
  tabpanel_ns <- NS(namespace=tabpanel_name)
  tabPanel(title='Score in clusters', value='scores_in_cluster',
    column(width=8, 
           div(style="display: inline-block;vertical-align:top; width: 150px;", cluster_id_picker.ui(id=tabpanel_name, opts=list(multiple=FALSE))),
           div(style="display: inline-block;vertical-align:top; width: 150px;", feature_picker.ui(id=tabpanel_name, gene_modules_opts=list(multiple=TRUE), include_feature_type=FALSE, include_values_range=FALSE)),
           plotOutput(outputId=tabpanel_ns('scores_in_cluster_ridges')) %>% withSpinner()),
    column(width=4, boxPlus(title='Options', width=12, closable=FALSE,
                            cluster_resolution_picker.ui(id=tabpanel_name)),
                            # add_gene_module.ui(id=tab)),
                    boxPlus(title='Map overview', width=12, closable=FALSE,
                            reduced_dimension_plot.ui(id=tabpanel_name, feature='selected_cluster_ids'),
                            reduction_method_picker.ui(id=tabpanel_name)),
                    boxPlus(title='Gene modules', width=12, closable=FALSE,
                            show_genes_in_module.ui(id=tabpanel_name)))) -> scores_in_cluster
  get0(env=module_servers_to_call, x=tabpanel_name) %>% append(sprintf(fmt='%s.server', module)) %>% assign(env=module_servers_to_call, x=tabpanel_name)

  tabBox(id=NS(namespace=id, id='gene_modules_scores_plot_box'),
         width=12, height='1250px',
         score_in_clusters, scores_in_cluster) -> ridges_plot_box

  # return ui element(s)
  tagList(ridges_plot_box)
}

#' Produce the ggplot object for a gene module score
#' 
#' @import ggridges
#' 
#' @rdname gene_module_score_plot
#'
gene_module_score_plot.server <- function(input, output, session, seurat, ...) {
  session$ns('') %>% sprintf(fmt='### %sgene_module_score_plot.server') %>% message()

  # get environments containing variables to run/configure this object
  collect_environments(id=parent.frame()$id, module='gene_module_score_plot') # provides `seuratvis_env`, `server_env` and `module_env`
  tab <- parent.frame()$id
  session_server <- get(x='session', env=server_env)
  input_server <- get(x='input', env=server_env)

  # restrict the feature selection to gene modules only
  updateRadioGroupButtons(session=session_server, inputId=NS(namespace=str_c(tab, 'score_in_clusters', sep='_'), id='feature_type'),
                          choices=list(`Gene modules`='gene_modules'),
                          selected='gene_modules')
  updateRadioGroupButtons(session=session_server, inputId=NS(namespace=str_c(tab, 'scores_in_cluster', sep='_'), id='feature_type'),
                          choices=list(`Gene modules`='gene_modules'),
                          selected='gene_modules')

  # react to the tab panel switching
  observeEvent(eventExpr=input$gene_modules_scores_plot_box, ignoreInit=FALSE, handlerExpr={})

  # render the gene module score in multiple clusters ggridges plot
  renderPlot(expr={
    tabpanel <- str_c(tab, input$gene_modules_scores_plot_box, sep='_')

    # make sure these elements are defined
    req(seurat$picked_cluster_resolution_idents[[tabpanel]])
    req(seurat$picked_feature_values[[tabpanel]])

    # send a message
    session$ns('') %>% sprintf(fmt='### %sgene_module_score_plot.server-renderPlot-IDs') %>% message()
   
    # make variables for shorthand    
    picked_gene_module <- input$picked_feature
    selected_cluster_idents <- input$cluster_id_picker

    # plot cluster (y) and score(x)
    cbind(seurat$picked_cluster_resolution_idents[[tabpanel]],
          seurat$picked_feature_values[[tabpanel]]) %>%
      filter(ident %in% selected_cluster_idents) %>%
      set_names(c('cluster_id', 'module_score')) %>%
      ggplot() +
      aes(x=module_score, y=cluster_id, fill=cluster_id) +
      labs(x='Gene module score', y='Cluster ID', fill='Cluster ID') +
      geom_density_ridges(alpha=0.8, size=0.7, colour='black') +
      # scale_fill_cyclical(values=c('grey50','grey70')) +
      theme_bw() +
      theme(axis.text.y=element_text(vjust=0),
            legend.position='none',
            panel.grid.major.x=element_blank(),
            panel.grid.minor.x=element_blank())}, height=1100) -> output$score_in_clusters_ridges

  # render the gene modules score in a cluster ggridges plot
  renderPlot(expr={
    tabpanel <- str_c(tab, input$gene_modules_scores_plot_box, sep='_')

    # make sure these elements are defined
    req(seurat$picked_cluster_resolution_idents[[tabpanel]])
    req(seurat$picked_feature_values[[tabpanel]])

    # send a message
    session$ns('') %>% sprintf(fmt='### %sgene_module_score_plot.server-renderPlot-GMs') %>% message()
   
    # make variables for shorthand    
    picked_gene_modules <- input$picked_feature
    selected_cluster_ident <- input$cluster_id_picker

    # plot cluster (y) and score(x)
    cbind(seurat$picked_cluster_resolution_idents[[tabpanel]],
          seurat$picked_feature_values[[tabpanel]]) %>%
      filter(ident==selected_cluster_ident) %>%
      gather(key=gene_module, value=module_score, -ident) %>%
      ggplot() +
      aes(x=module_score, y=gene_module) +
      labs(x='Gene module score', y='Gene module') +
      geom_density_ridges(alpha=0.8, size=0.7, colour='black') +
      scale_fill_cyclical(values=c('grey50','grey70')) +
      theme_bw() +
      theme(axis.text.y=element_text(vjust=0),
            legend.position='none',
            panel.grid.major.x=element_blank(),
            panel.grid.minor.x=element_blank())}, height=1100) -> output$scores_in_cluster_ridges


}
