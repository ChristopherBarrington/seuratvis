#'
#' 
dimensionality.plot <- function(id)
  plotOutput(outputId=NS(id, 'plot')) %>% withSpinner()

#'
#' 
dimensionality.jackstraw_pvalue <- function(input, output, session, seurat, picked_reduction)
  renderPlot(expr={
    # make sure these elements are defined
    req(seurat$object)
    req(picked_reduction$method)

    # make a plot
    as.data.frame(Seurat::JS(object=seurat$object[[picked_reduction$method]], slot='overall')) %>%
      ggplot() +
      aes(x=PC, y=-log10(Score)) +
      labs(x='Principle component', y='-log10(score)') +
      geom_smooth(method='loess') +
      geom_point(shape=4) +
      theme_bw() +
      theme(legend.background=element_blank(),
            legend.justification=c(1,1),
            legend.position=c(1,1))}) -> output$plot

#'
#' 
dimensionality.jackstraw <- function(input, output, session, seurat, picked_reduction, picked_components)
  renderPlot(expr={
    req(seurat$object)
    req(picked_reduction$method)
    req(picked_components$picked)

    # make variables for shorthand    
    object <- seurat$object
    reduction_name <- picked_reduction$method
    components_range <- picked_components$picked
    min_component <- min(components_range)
    max_component <- max(components_range)

    Seurat::JS(object=object[[reduction_name]], slot='empirical') %>%
      as.data.frame() %>%
      rownames_to_column('Contig') %>%
      gather(key='PC', value='Value', -Contig) %>%
      mutate(PC={stringr::str_remove(PC, '^PC') %>% as.numeric()}) %>%
      filter(dplyr::between(PC, left=min_component, right=max_component)) %>%
      left_join(y=as.data.frame(Seurat::JS(object=object[[reduction_name]], slot='overall')), by='PC') %>%
      mutate(PC_colour=sprintf('PC %d: %1.3g', PC, Score)) %>%
      mutate(PC_colour=factor(PC_colour, levels={unique(PC_colour) %>% str_sort(numeric=TRUE)})) %>%
      ggplot() +
      aes(sample=Value, colour=PC_colour) +
      labs(x='Theoretical [runif(1000)]', y='Empirical', colour='PC pvalue') +
      stat_qq(distribution=qunif, alpha=0.5) +
      geom_abline(intercept=0, slope=1, linetype='dashed') +
      coord_flip() +
      theme_bw() +
      theme(legend.position='none')}) -> output$plot

#'
#' 
dimensionality.elbow <- function(input, output, session, seurat, picked_reduction)
  renderPlot(expr={
    req(seurat$object)
    req(picked_reduction$method)

    data.frame(Y=Stdev(object=seurat$object, reduction=picked_reduction$method)) %>%
      mutate(X=seq(n())) -> data

    stdev <- Stdev(object=seurat$object, reduction=picked_reduction$method)
    pct <- stdev / sum(stdev) * 100 # Determine percent of variation associated with each PC
    cumu <- cumsum(pct) # Calculate cumulative percents for each PC
    
    co1 <- which(cumu > 90 & pct < 5)[1] # Determine which PC exhibits cumulative percent greater than 90% and % variation associated with the PC as less than 5
    co2 <- sort(which((pct[1:length(pct) - 1] - pct[2:length(pct)]) > 0.1), decreasing=TRUE)[1] + 1 # Determine the difference between variation of PC and subsequent PC

    data.frame(co1=co1, co2=co2) %>%
      mutate(pcs=pmin(co1, co2),
             min=pmin(co1, co2),
             max=pmax(co1, co2)) -> dimensions

    data %>%
      mutate(selected=cut(x=X, breaks=c(0, co1, co2, length(stdev)+1))) %>%
      ggplot(data=.) +
      aes(x=X, y=Y, colour=selected) +
      labs(x='Principal component', y='Standard deviation', colour='Selected') +
      geom_line(colour='grey', size=1) +
      geom_point(size=2) +
      theme_bw() +
      theme(aspect.ratio=1,
            legend.position='bottom',
            legend.title=element_blank(),
            panel.grid.minor=element_blank(),
            strip.background=element_rect(fill=NA))}) -> output$plot
#'
#' 
dimensionality.top_features_pca_heatmap <- function(input, output, session, seurat, picked_reduction, picked_components)
  renderPlot(expr={
    # make sure these elements are defined
    req(seurat$object)
    req(picked_reduction$method)
    req(picked_components$picked)
   
    # make variables for shorthand    
    object <- seurat$object
    reduction_name <- picked_reduction$method
    selected_component <- picked_components$picked

    # render the heatmap
    #! TODO: make this a nice ggplot
    if(!DefaultAssay(object=object[[reduction_name]]) %in% Assays(object)) {
      ggplot()+aes()+annotate(geom='text', label='Nothing to see here', x=0, y=0)+theme_void()
    } else {
      DefaultAssay(object=object) <- DefaultAssay(object=object[[reduction_name]])
      DimHeatmap(object=object, reduction=reduction_name,
                 dims=as.numeric(selected_component),
                 disp.min=-2.5, disp.max=2.5,
                 cells=2000, balanced=TRUE, fast=FALSE)}}) -> output$plot
