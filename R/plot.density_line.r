#' 
#' 
density.plot <- function(id, feature) {
  plotOutput(outputId=NS(id, 'density'), brush=brushOpts(id=NS(id, 'brush'), direction='x')) %>% withSpinner()
}

#'
#'
density_plot.add_common <- function(ggplot_object)
  ggplot_object +
    theme_bw() +
    theme()

#'
#'
density_plot.n_features.server <- function(input, output, session, seurat) {
  plot_brush <- reactiveValues()

  renderPlot(expr={
    req(seurat$n_features_values)

    data.frame(y=seurat$n_features_values) %>%
      {ggplot(data=.) +
       aes(x=y) +
       labs(x='Features detected per cell', y='Density') +
       stat_density(geom='line', trim=TRUE, size=2) +
       scale_x_log10(breaks=major_breaks_log10, minor_breaks=minor_breaks_log10, labels=function(x) scales::comma(x, accuracy=1))} %>%
      density_plot.add_common()}) -> output$density

  # update plot_brush when the brush is used
  observe(label='density_plot/n_features/brush', x={
    req(input$brush)

    # get brush positions, round up or down
    low <- floor(input$brush$xmin)
    high <- ceiling(input$brush$xmax)

    # update the reactive
    plot_brush$min <- low
    plot_brush$max <- high})

  # return the reactiveValues
  return(plot_brush)
}

#'
#'
density_plot.n_umi.server <- function(input, output, session, seurat) {
  plot_brush <- reactiveValues()

  renderPlot(expr={
    req(seurat$n_umi_values)

    data.frame(y=seurat$n_umi_values) %>%
      {ggplot(data=.) +
       aes(x=y) +
       labs(x='Total UMI per cell', y='Density') +
       stat_density(geom='line', trim=TRUE, size=2) +
       scale_x_log10(breaks=major_breaks_log10, minor_breaks=minor_breaks_log10, labels=function(x) scales::comma(x, accuracy=1))} %>%
      density_plot.add_common()}) -> output$density

  # update plot_brush when the brush is used
  observe(label='density_plot/n_umi/brush', x={
    req(input$brush)

    # get brush positions, round up or down
    low <- floor(input$brush$xmin)
    high <- ceiling(input$brush$xmax)

    # update the reactive
    plot_brush$min <- low
    plot_brush$max <- high})

  # return the reactiveValues
  return(plot_brush)
}

#'
#'
density_plot.proportion_mt.server <- function(input, output, session, seurat) {
  plot_brush <- reactiveValues()

  renderPlot(expr={
    req(seurat$proportion_mt_values)

    data.frame(y=seurat$proportion_mt_values) %>%
      {ggplot(data=.) +
       aes(x=y) +
       labs(x='Proportion mitochondrial expression', y='Density') +
       stat_density(geom='line', trim=TRUE, size=2)} %>%
      density_plot.add_common()}) -> output$density

  # update plot_brush when the brush is used
  observe(label='density_plot/proportion_mt/brush', x={
    req(input$brush)

    # get brush positions, round up or down
    high <- round(input$brush$xmax+0.05, digits=1)

    # update the reactive
    plot_brush$max <- high})

  # return the reactiveValues
  return(plot_brush)
}







###  # react to the brush
###  #! TODO: these should not need session servers?
###  observeEvent(eventExpr=input$brush, handlerExpr={
###    # send a message
###    session$ns('') %>% sprintf(fmt='### %sdensity_plot.server-observeEvent-input$brush [%s]', input$brush) %>% message()
###
###    # update feature-specific reactives and ui elements
###    if(module_env$feature=='nFeature_RNA') {
###      # get brush positions, round up or down
###      low <- floor(input$brush$xmin)
###      high <- ceiling(input$brush$xmax)
###
###      # update the ui
###      for(id in module_environments$features_per_cell_filters$id) {
###        updateNumericInput(session=session_server, inputId=NS(id, 'min_features'), value=low)
###        updateNumericInput(session=session_server, inputId=NS(id, 'max_features'), value=high)
###      }
###
###      # update the reactive
###      cell_filtering$features_per_cell_min <- low
###      cell_filtering$features_per_cell_max <- high
###    } else if(module_env$feature=='nCount_RNA') {
###      # get brush positions, round up or down
###      low <- floor(input$brush$xmin)
###      high <- ceiling(input$brush$xmax)
###
###      # update the ui
###      for(id in module_environments$total_umi_per_cell_filters$id) {
###        updateNumericInput(session=session_server, inputId=NS(namespace=id, id='min_umis'), value=low)
###        updateNumericInput(session=session_server, inputId=NS(namespace=id, id='max_umis'), value=high)
###      }
###
###      # update the reactive
###      cell_filtering$total_umi_per_cell_min <- low
###      cell_filtering$total_umi_per_cell_max <- high
###    } else if(module_env$feature=='percent_mt') {
###      # get the max position of the brush to 1dp
###      high <- round(input$brush$xmax+0.05, digits=1)
###
###      # update the ui
###      #! TODO remove this dependency, save the value to a reactive anf react to it in the gui module
###      for(id in module_environments$percent_mt_per_cell_filters$id)
###        updateNumericInput(session=session_server, inputId=NS(namespace=id, id='max_percent_mt'), value=high)
###
###      # update the reactive
###      cell_filtering$max_percent_mitochondria <- high
###    }})
###
###  # reset the brush when n features variable is changed
###  #! TODO: these should not require the NS() call ...
###  observeEvent(eventExpr=seurat$reset_n_features, handlerExpr={
###    # send a message
###    session$ns('') %>% sprintf(fmt='### %sdensity_plot.server-observeEvent-seurat$reset_n_features [%s]', seurat$reset_n_features) %>% message()
###    
###    session$resetBrush(NS(namespace=id, id= 'brush'))})
###
###  # reset the brush when n umi variable is changed
###  observeEvent(eventExpr=seurat$reset_n_umi, handlerExpr={
###    # send a message
###    session$ns('') %>% sprintf(fmt='### %sdensity_plot.server-observeEvent-seurat$reset_n_umi [%s]', seurat$reset_n_umi) %>% message()
###    
###    session$resetBrush(NS(namespace=id, id='brush'))})
###
###  # reset the brush when proportion mitochondrial UMI variable is changed
###  observeEvent(eventExpr=seurat$reset_proportion_mt, handlerExpr={
###    # send a message
###    session$ns('') %>% sprintf(fmt='### %sdensity_plot.server-observeEvent-seurat$reset_proportion_mt [%s]', seurat$reset_proportion_mt) %>% message()
###    
###    session$resetBrush(NS(namespace=id, id='brush'))})
###}
