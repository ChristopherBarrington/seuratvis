#'
#'
feature_ridges.plot <- function(id)
  plotOutput(outputId=NS(id, 'ridges'), height='70vh') %>% withSpinner()

#' Plot ridges for a feature across all idents
#' 
#' @import ggridges
#' 
feature_ridge_by_idents.server <- function(input, output, session, picked_feature, picked_clusters) {
  renderPlot(expr={
    if(!isTruthy(picked_clusters$idents) | !isTruthy(picked_feature$values))
      return(missing_data_plot())

    cbind(ident=picked_clusters$idents,
          value=picked_feature$values) %>%
      filter(ident %in% picked_clusters$picked_idents) %>%
      mutate(ident=factor(ident, levels=picked_clusters$distinct_idents)) %>%
      set_names(c('cluster_id', 'module_score')) %>%
      ggplot() +
      aes(x=module_score, y=cluster_id, fill=cluster_id) +
      labs(x='Gene module score', y='Cluster ID', fill='Cluster ID') +
      geom_density_ridges(alpha=0.8, size=0.7, colour='black') +
      scale_fill_discrete(drop=FALSE) +
      theme_bw() +
      theme(axis.text.y=element_text(vjust=0),
            legend.position='none',
            panel.grid.major.x=element_blank(),
            panel.grid.minor.x=element_blank())}) -> output$ridges
}

#' Plot ridges for multiple features in an ident class
#' 
#' @import ggridges
#' 
feature_ridges_by_ident.server <- function(input, output, session, picked_feature, picked_clusters) {
  renderPlot({
    if(!isTruthy(picked_clusters$idents) | !isTruthy(picked_clusters$picked_idents) | !isTruthy(picked_feature$values))
      return(missing_data_plot())

    cbind(ident=picked_clusters$idents,
          picked_feature$values) %>%
      filter(ident==picked_clusters$picked_idents) %>%
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
            panel.grid.minor.x=element_blank())}) -> output$ridges
}