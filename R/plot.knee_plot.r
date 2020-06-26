#' 
#' 
knee.plot <- function(id) {
  plotOutput(outputId=NS(id, 'knee')) %>% withSpinner()
}

#'
#' 
knee_plot.add_common <- function(ggplot_object, max_value)
  ggplot_object +
    labs(x='Ranked cells', colour='Cells passing threshold') +
    geom_hline(yintercept=max_value, size=1) +
    geom_point(alpha=1, shape=16) +
    scale_x_log10(breaks=major_breaks_log10, minor_breaks=minor_breaks_log10, labels=function(x) scales::comma(x, accuracy=1)) +
    scale_colour_brewer(palette='Set1', direction=1) +
    theme_bw() +
    theme(legend.position='none', legend.title=element_blank())

#' 
#'
knee_plot.n_features.server <- function(input, output, session, seurat, cell_filtering) {
  renderPlot(expr={
    req(cell_filtering$n_features_min)
    req(cell_filtering$n_features_max)
    req(seurat$n_features_values)

    min_value <- cell_filtering$n_features_min
    max_value <- cell_filtering$n_features_max
    
    # start the knee plot
    data.frame(y=seurat$n_features_values) %>%
      arrange(desc(y)) %>%
      mutate(x=seq(n()),
             pass=between(x=y, left=min_value, right=max_value)) %>%
      {ggplot(data=.) +
       aes(x=x, y=y, colour=pass) +
       labs(y='Features detected per cell') +
       geom_hline(yintercept=min_value, size=1) +
       scale_y_log10(breaks=major_breaks_log10, minor_breaks=minor_breaks_log10, labels=function(y) scales::comma(y, accuracy=1))} %>%
      knee_plot.add_common(max_value=max_value)}) -> output$knee
}

#' 
#'
knee_plot.n_umi.server <- function(input, output, session, seurat, cell_filtering) {
  renderPlot(expr={
    req(cell_filtering$n_umi_min)
    req(cell_filtering$n_umi_max)
    req(seurat$n_umi_values)

    min_value <- cell_filtering$n_umi_min
    max_value <- cell_filtering$n_umi_max
    
    # start the knee plot
    data.frame(y=seurat$n_umi_values) %>%
      arrange(desc(y)) %>%
      mutate(x=seq(n()),
             pass=between(x=y, left=min_value, right=max_value)) %>%
      {ggplot(data=.) +
       aes(x=x, y=y, colour=pass) +
       labs(y='Total UMI per cell') +
       geom_hline(yintercept=min_value, size=1) +
       scale_y_log10(breaks=major_breaks_log10, minor_breaks=minor_breaks_log10, labels=function(y) scales::comma(y, accuracy=1))} %>%
      knee_plot.add_common(max_value=max_value)}) -> output$knee
}

#' 
#'
knee_plot.proportion_mt.server <- function(input, output, session, seurat, cell_filtering) {
  renderPlot(expr={
    req(cell_filtering$proportion_mt_max)
    req(seurat$proportion_mt_values)

    max_value <- cell_filtering$proportion_mt_max

    # start the knee plot
    data.frame(y=seurat$proportion_mt_values) %>%
      arrange(desc(y)) %>%
      mutate(x=seq(n()),
             pass=y<=max_value) %>%
      {ggplot(data=.)+
       aes(x=x, y=y, colour=pass)+
       labs(y='Proportion mitochondrial expression')} %>%
      knee_plot.add_common(max_value=max_value)}) -> output$knee
}
