#'
#' 
dimension_reduction.plot <- function(id) 
  tagList(NS(id, 'map') %>% plotOutput() %>% withSpinner())

#'
#' 
dimension_reduction.plotbox <- function(id, n=0) {
  # make the ui elements
  ## drop down config
  dropdownButton(inputId=NS(id,'ddn'),
                 circle=TRUE, status='primary', icon=icon('wrench'), size='xs',
                 uiOutput(outputId=NS(id,'dropdown'))) -> dropdown_menu

  tags$head(tags$style(HTML(sprintf('#%s {background: #222D32 !important;}', NS(NS('dropdown-menu',id),'ddn'))))) -> css

  ## map plot output
  NS(id, 'map') %>% plotOutput() %>% withSpinner() -> plot_output

  # return the box
  tagList(div(id=str_c('boxid', n), boxPlus(title=textOutput(outputId=NS(id, 'box_title')), status='primary', closable=FALSE, width=4, dropdown_menu, css, plot_output)))
}

#'
#' 
dimension_reduction.show_cluster_idents.server <- function(input, output, session, dimension_reduction, picked_colours, opacity, point_size, cluster_resolution) {
  renderPlot(bg='transparent', expr={
    req(dimension_reduction$embeddings)
    req(cluster_resolution$idents)
 
    # prepare the data
    cbind(dimension_reduction$embeddings %>% set_names(c('X','Y')),
          ident=cluster_resolution$idents) %>%
      sample_frac(size=1) -> data
 
    # make a ggplot
    ggplot(data=data) +
      aes(x=X, y=Y, colour=ident) +
      geom_hline(yintercept=0, colour='grey90') + geom_vline(xintercept=0, colour='grey90') +
      geom_point(size=point_size$size, alpha=opacity$alpha) +
      theme_void() +
      theme(aspect.ratio=1,
            legend.position='none',
            panel.background=element_rect(fill=picked_colours$background)) -> map
 
    # if labels should be included, add them here
    if(cluster_resolution$label_clusters) {
      # find the mean positions of the clusters
      data %>%
        group_by(ident) %>%
        summarise(X=mean(X), Y=mean(Y)) -> data_labels
 
      # add labels to the plot
      map +
        geom_label_repel(data=data_labels, mapping=aes(label=ident),
                         colour='black', size=12/(14/5)) -> map
    }
 
    # return the map
    map}) -> output$map
}

#'
#'
dimension_reduction.show_selected_clusters.server <- function(input, output, session, dimension_reduction, picked_colours, point_size, cluster_resolution) {
  renderPlot(bg='transparent', expr={
    req(dimension_reduction$embeddings)
    req(cluster_resolution$idents)
 
    # prepare the data and make the plot
    cbind(dimension_reduction$embeddings %>% set_names(c('X','Y')),
          ident=cluster_resolution$idents) %>%
      mutate(is_selected=ident %in% cluster_resolution$picked_idents) %>%
      arrange(is_selected) %>%
      ggplot() +
      aes(x=X, y=Y, colour=ident, alpha=is_selected) +
      geom_hline(yintercept=0, colour='grey90') + geom_vline(xintercept=0, colour='grey90') +
      geom_point(size=point_size$size) +
      scale_alpha_manual(values=c(`FALSE`=0.05, `TRUE`=1)) +
      theme_void() +
      theme(aspect.ratio=1,
            legend.position='none',
            panel.background=element_rect(fill=picked_colours$background)) -> map

    # return the plot
    map}) -> output$map
}

#'
#' 
dimension_reduction.highlight_feature.server <- function(input, output, session, dimension_reduction, picked_feature, picked_colours, opacity, point_size) {
  # render the box title if defined
  renderText({ifelse(is.null(picked_feature$name), 'Select feature using dropdown', picked_feature$name)}) -> output$box_title

  # render the plot
  renderPlot(bg='transparent', expr={
    req(dimension_reduction$embeddings)
    req(picked_feature$values)

    # prepare the data and start the plot
    cbind(dimension_reduction$embeddings %>% set_names(c('X','Y')),
          picked_feature$values %>% set_names('value')) %>%
      arrange(value) %>%
      ggplot() +
      aes(x=X, y=Y, colour=value) +
      geom_hline(yintercept=0, colour='grey90') + geom_vline(xintercept=0, colour='grey90') +
      geom_point(size=point_size$size, alpha=opacity$alpha) +
      theme_void() +
      theme(aspect.ratio=1,
            legend.position='none',
            panel.background=element_rect(fill=picked_colours$background, colour='black'),
            strip.background=element_rect(fill=picked_colours$background, colour='black')) -> map

    # if the feature is numeric, colour the points otherwise facet the plot
    if(is.numeric(map$data$value)) {
      # get the colour values and range
      c_low <- picked_colours$low
      c_mid <- picked_colours$mid
      c_high <- picked_colours$high
      c_range <- picked_feature$values_range

      # make a colour scale
      colour_gradient <- scale_colour_gradient(low=c_low, high=c_high, limits=c_range, oob=scales::squish)

      # if the values cross zero, make a new colour scale
      if(c_range %>% sign() %>% Reduce(f='*') %>% magrittr::equals(-1)) {
         c_low <- 'blue'
         c_high <- 'red'
         colour_gradient <- scale_colour_gradientn(colours=c(low=c_low, mid=c_mid, high=c_high), 
                                                  values={c_range %>% c(0) %>% sort() %>% scales::rescale()},
                                                  limits=c_range, breaks=0)
      }

      # add the colour scale and a legend
      map +
        colour_gradient +
        # labs(colour=picked_feature$name) +
        guides(colour=guide_colourbar(title.position='top',
                                      title.hjust=1,
                                      frame.colour='black',
                                      ticks=TRUE)) +
        theme(legend.direction='horizontal',
              legend.justification=c(1,0),
              legend.key.height=unit(0.5,'line'),
              legend.key.width=unit(0.5,'line'),
              legend.position=c(1,0),
              legend.text=element_blank(),
              legend.title=element_blank(),
              legend.title.align=0,
              panel.background=element_rect(fill=picked_colours$background, colour='black')) -> map
    } else {
      # facet the map by the factor levels
      map <- map + facet_wrap(~value, scales='fixed')
    }

    # return the plot
    map}) -> output$map
}

#'
#' 
dimension_reduction.split_highlight_feature.server <- function(input, output, session, dimension_reduction, picked_feature, picked_colours, opacity, point_size, picked_splitter) {
  # render the plot
  renderPlot(bg='transparent', expr={
    req(dimension_reduction$embeddings)
    req(picked_feature$values)
    req(picked_splitter$values)

    # prepare the data and start the plot
    cbind(dimension_reduction$embeddings %>% set_names(c('X','Y')),
          picked_feature$values %>% set_names('value'),
          picked_splitter$values %>% set_names('splitter')) %>%
      (function(x) rbind(x, mutate(x, splitter='unsplit dataset'))) %>%
      arrange(value) %>%
      ggplot() +
      aes(x=X, y=Y, colour=value) +
      geom_hline(yintercept=0, colour='grey90') + geom_vline(xintercept=0, colour='grey90') +
      geom_point(size=point_size$size, alpha=opacity$alpha) +
      facet_wrap(~splitter) +
      theme_void() +
      theme(aspect.ratio=1,
            legend.position='none',
            panel.background=element_rect(fill=picked_colours$background, colour='black'),
            strip.background=element_rect(fill='white', colour='black'),
            strip.text=element_text(size=15)) -> map

    # if the feature is numeric, colour the points otherwise facet the plot
    if(is.numeric(map$data$value)) {
      # get the colour values and range
      c_low <- picked_colours$low
      c_mid <- picked_colours$mid
      c_high <- picked_colours$high
      c_range <- picked_feature$values_range

      # make a colour scale
      colour_gradient <- scale_colour_gradient(low=c_low, high=c_high, limits=c_range, oob=scales::squish)

      # if the values cross zero, make a new colour scale
      if(c_range %>% sign() %>% Reduce(f='*') %>% magrittr::equals(-1)) {
         c_low <- 'blue'
         c_high <- 'red'
         colour_gradient <- scale_colour_gradientn(colours=c(low=c_low, mid=c_mid, high=c_high), 
                                                  values={c_range %>% c(0) %>% sort() %>% scales::rescale()},
                                                  limits=c_range, breaks=0)
      }

      # add the colour scale and a legend
      map +
        colour_gradient +
        # labs(colour=picked_feature$name) +
        guides(colour=guide_colourbar(title.position='top',
                                      title.hjust=1,
                                      frame.colour='black',
                                      ticks=TRUE)) +
        theme(legend.direction='horizontal',
              legend.justification=c(1,0),
              legend.key.height=unit(0.5,'line'),
              legend.key.width=unit(0.5,'line'),
              legend.position=c(1,0),
              legend.text=element_blank(),
              legend.title=element_blank(),
              legend.title.align=0,
              panel.background=element_rect(fill=picked_colours$background, colour='black'),
              strip.background=element_rect(fill='white', colour='black')) -> map
    } else {
      # facet the map by the factor levels
      map <- map + facet_wrap(~value, scales='fixed')
    }

    # return the plot
    map}, height=900) -> output$map
}
