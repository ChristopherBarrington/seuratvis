#' 
#' 
feature_value_per_cluster.plot <- function(id)
  plotOutput(outputId=NS(id,'boxplot')) %>% withSpinner()

#' 
#'
feature_value_per_cluster.server <- function(input, output, session, picked_feature, cluster_resolution, picked_colours) {
  renderPlot(bg='transparent', expr={
    req(cluster_resolution$idents)
    req(picked_feature$values)

    # get the data to plot
    cbind(ident=cluster_resolution$idents,
          picked_feature$values) %>%
      mutate(x={as.character(ident) %>% as.numeric()}) -> data

    # if the picked feature has numeric values
    if(is.numeric(data$value)) {
      # make a values range-type plot
      data %>%
        # filter(value>0) %>%
        group_by(ident, x) %>%
        summarise(q25=quantile(value, 0.25), q75=quantile(value, 0.75), median=median(value)) %>%
        mutate(is_selected_cluster_id=ident %in% input$cluster_id_picker) %>%
        mutate(iqr=q75-q25, lower=q25-1.5*iqr, upper=q75+1.5*iqr) -> cluster_data_summary

      cluster_data_summary %>%
        gather(key=key, value=y, lower, upper) %>%
        ggplot() +
        aes(x=x, y=y, colour=ident) +
        labs(y='Feature value (median ± 1.5x IQR)') +
        geom_line(size=1) +
        geom_point(mapping=aes(y=median, fill=is_selected_cluster_id), shape=21, size=3, colour='black') +
        scale_fill_manual(values=c(`FALSE`='grey90', `TRUE`='darkorange1')) -> feature_plot
    } else { # the picked feature is factor-like
      # make a frequency bar plot
      data %>%
        ggplot() +
        aes(x=x, fill=value) +
        labs(y='Frequency') +
        geom_bar(position='dodge') -> feature_plot
    }

    # add the shared plotting elements
    feature_plot +
      labs(x='Cluster identifier') +
      scale_x_continuous(breaks=seq(from=0, to=100, by=2), minor_breaks=seq(from=1, to=99, by=2), sec.axis=dup_axis(name=NULL)) +
      theme_bw() +
      theme(legend.position='none',
            panel.background=element_rect(fill=picked_colours$background, colour='black'),
            strip.background=element_rect(fill=picked_colours$background, colour='black'))}) -> output$boxplot
}
