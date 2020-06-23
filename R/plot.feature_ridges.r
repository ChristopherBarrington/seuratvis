#'
#'
feature_ridges.plot <- function(id)
  plotOutput(outputId=NS(id, 'ridges'), height='70vh') %>% withSpinner()

#'
#' @import ggridges
#' 
feature_ridge_by_idents.server <- function(input, output, session, picked_feature, picked_clusters) {
  renderPlot(expr={
    # req(picked_clusters$idents)
    # req(picked_feature$values)
    if(!isTruthy(picked_clusters$idents) | !isTruthy(picked_feature$values))
      return(NULL)

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

#'
#' 
feature_ridges_by_ident <- function() {}
