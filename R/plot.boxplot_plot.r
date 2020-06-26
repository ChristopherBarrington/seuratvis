#' 
#' 
boxplot.plot <- function(id) {
  plotOutput(outputId=NS(id, 'boxplot')) %>% withSpinner()
}

#'
#' 
boxplot_plot.add_common <- function(ggplot_object)
  ggplot_object +
    geom_point(shape=16, size=0.6, colour='lightgrey', alpha=0.6, position=position_jitter(width=0.5)) +
    geom_boxplot(fill=NA, width=0.4, size=0.9, outlier.size=0.6) +
    coord_cartesian(xlim=c(-1, 1)) +
    theme_bw() +
    theme(axis.ticks.length.x=unit(0, 'pt'),
          axis.text.x=element_text(colour='white'),
          axis.title.x=element_text(colour='white'),
          panel.grid.major.x=element_blank(),
          panel.grid.minor.x=element_blank())

#' 
#'
boxplot_plot.n_features.server <- function(input, output, session, seurat, cell_filtering) {
  renderPlot(expr={
    req(cell_filtering$n_features_min)
    req(cell_filtering$n_features_max)
    req(seurat$n_features_values)

    # get thresholds of detected features
    min_y <- cell_filtering$n_features_min
    max_y <- cell_filtering$n_features_max
    
    # start the boxplot
    data.frame(y=seurat$n_features_values) %>%
      {ggplot(data=.) +
       aes(x=0, y=y) +
       labs(y='Features detected per cell') +
       annotate(geom='rect', xmin=-Inf, xmax=Inf, ymin=min_y, ymax=max_y, fill='orange', alpha=0.2) +
       scale_y_continuous(labels=function(y) scales::comma(y, accuracy=1))} %>%
      boxplot_plot.add_common()}) -> output$boxplot
}

#' 
#'
boxplot_plot.n_umi.server <- function(input, output, session, seurat, cell_filtering) {
  renderPlot(expr={
    req(cell_filtering$n_umi_min)
    req(cell_filtering$n_umi_max)
    req(seurat$n_umi_values)

    # get thresholds of total UMI
    min_y <- cell_filtering$n_umi_min
    max_y <- cell_filtering$n_umi_max
    
    # start the boxplot
    data.frame(y=seurat$n_umi_values) %>%
      {ggplot(data=.) +
       aes(x=0, y=y) +
       labs(y='Total UMI per cell') +
       annotate(geom='rect', xmin=-Inf, xmax=Inf, ymin=min_y, ymax=max_y, fill='orange', alpha=0.2) +
       scale_y_continuous(labels=function(y) scales::comma(y, accuracy=1))} %>%
     boxplot_plot.add_common()}) -> output$boxplot
}

#' 
#'
boxplot_plot.proportion_mt.server <- function(input, output, session, seurat, cell_filtering) {
  renderPlot(expr={
    req(cell_filtering$proportion_mt_max)
    req(seurat$proportion_mt_values)

    # get max fraction Mt UMI
    max_y <- cell_filtering$proportion_mt_max

    # start the boxplot
    data.frame(y=seurat$proportion_mt_values) %>%
      {ggplot(data=.) +
       aes(x=0, y=y) +
       labs(y='Proportion mitochondrial expression') +
       annotate(geom='rect', xmin=-Inf, xmax=Inf, ymin=0, ymax=max_y, fill='orange', alpha=0.2)} %>%
      boxplot_plot.add_common()}) -> output$boxplot
}
