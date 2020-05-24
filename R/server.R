
10^(0:9) -> major_breaks_log10
(2:9) * rep(major_breaks_log10, each=8) -> minor_breaks_log10

module_environments <- new.env()
seurat_object.reactions <- reactiveValues()
filtering_parameters.reactions <- reactiveValues()

shinyAppServer <- function(input, output, session) {

  progress <- shiny::Progress$new(session=session, min=0, max=4/10)
  on.exit(progress$close())
  progress$set(value=0, message='Loading environment')

  # ###############################################################################################
  # scour the session for Seurat objects and populate the UI --------------------------------------
  available_seurat_objects <- find_seurat_objects() # has to be here to access the global environment ... ?
  callModule(module=available_seurats.server, id='load_dataset')
  callModule(module=load_a_seurat.server, id='load_dataset')

  # ###############################################################################################
  # cell filtering tab ----------------------------------------------------------------------------

  ## react to cell filtering parameters
  reactiveValues(filtered_cell_set=NULL, n_cells=0, total_reads=0, median_reads_per_cell=0, median_genes_per_cell=0,
                 min_genes_per_cell=NULL, max_genes_per_cell=NULL,
                 min_expression_per_cell=NULL, max_expression_per_cell=NULL,
                 max_percent_mitochondria=NULL,
                 subset_conditions=list(),
                 subset_conditions.nCount='', subset_conditions.nFeature='', subset_conditions.percent_mt='') -> cell_filtering_data.reactions

  reactive(x={
    progress <- shiny::Progress$new(session=session, min=0, max=4/10)
    on.exit(progress$close())
    progress$set(value=0, message='Reacting to threshold parameters')

    progress$inc(detail='Extracting input options')
    min_expression_per_cell <- seurat_object.reactions$total_umi_per_cell_min
    max_expression_per_cell <- seurat_object.reactions$total_umi_per_cell_max
    min_genes_per_cell <- seurat_object.reactions$features_per_cell_min
    max_genes_per_cell <- seurat_object.reactions$features_per_cell_max
    max_percent_mitochondria <- seurat_object.reactions$percent_mt_per_cell_max

    progress$inc(detail='Filtering @meta.data')
    seurat_object.reactions$seurat@meta.data %>%
      filter(percent_mt <= max_percent_mitochondria &
             between(x=nFeature_RNA, left=min_genes_per_cell, right=max_genes_per_cell) &
             between(x=nCount_RNA, left=min_expression_per_cell, right=max_expression_per_cell)) -> filtered_cell_set

    progress$inc(detail='Saving results to reactiveValues')
    cell_filtering_data.reactions$filtered_cell_set <- filtered_cell_set
    cell_filtering_data.reactions$n_cells <- nrow(filtered_cell_set)
    cell_filtering_data.reactions$total_reads <- sum(filtered_cell_set$nCount_RNA)
    cell_filtering_data.reactions$median_reads_per_cell <- round(x=median(filtered_cell_set$nCount_RNA), digits=0)
    cell_filtering_data.reactions$median_genes_per_cell <- round(x=median(filtered_cell_set$nFeature_RNA), digits=0)

    # save values to filtering reactive
    filtering_parameters.reactions$n_cells <- nrow(filtered_cell_set)
    filtering_parameters.reactions$n_umi <- sum(filtered_cell_set$nCount_RNA)

    progress$inc(detail='Returning')
    NULL}) -> react_to_cell_filtering

  ## react to unique features detected density plot brush
  observeEvent(eventExpr=input$unique_genes_density.brush, handlerExpr={
    message('### input$unique_genes_density.brush')

    # update module elements
    low <- floor(input$unique_genes_density.brush$xmin)
    high <- ceiling(input$unique_genes_density.brush$xmax)

    for(id in module_environments$features_per_cell_filters$id) {
      updateTextInput(session=session, inputId=NS(id, 'min_features'), value=low)
      updateTextInput(session=session, inputId=NS(id, 'max_features'), value=high)
    }

    seurat_object.reactions$features_per_cell_min <- low
    seurat_object.reactions$features_per_cell_max <- high})

  ## react to total UMIs density plot brush
  observeEvent(eventExpr=input$total_expression_density.brush, handlerExpr={
    message('### input$total_expression_density.brush')

    # update module elements
    low <- floor(input$total_expression_density.brush$xmin)
    high <- ceiling(input$total_expression_density.brush$xmax)

    for(id in module_environments$total_umi_per_cell_filters$id) {
      updateTextInput(session=session, inputId=NS(namespace=id, id='min_umis'), value=low)
      updateTextInput(session=session, inputId=NS(namespace=id, id='max_umis'), value=high)
    }

    seurat_object.reactions$total_umi_per_cell_min <- low
    seurat_object.reactions$total_umi_per_cell_max <- high})

  ## react to total UMI per cell filters
  for(id in module_environments$total_umi_per_cell_filters$id)
    callModule(module=total_umi_per_cell_filter.server, id=id)

  ## react to features per cell filters
  for(id in module_environments$features_per_cell_filters$id)
    callModule(module=features_per_cell_filter.server, id=id)

  ## react to percent Mt per cell filters
  for(id in module_environments$percent_mt_per_cell_filters$id)
    callModule(module=percent_mt_per_cell_filter.server, id=id)

  ## react to show filtering parameters
  for(id in module_environments$show_filtering_parameters$id)
    callModule(module=show_filtering_parameters.server, id=id)

  ## react to percent mitochondria density plot brush
  observeEvent(eventExpr=input$percent_mitochondria_density.brush, handlerExpr={
    high <- round(input$percent_mitochondria_density.brush$xmax+0.05, digits=1)

    for(id in module_environments$percent_mt_per_cell_filters$id)
      updateTextInput(session=session, inputId=NS(namespace=id, id='max_percent_mt'), value=high)

    seurat_object.reactions$percent_mt_per_cell_max <- high})

  ## react to opening tab with a filtered object loaded
  observeEvent(input$sidebarmenu, {
    if(!is.null(seurat_object.reactions$seurat) & input$sidebarmenu=='cell_filtering-tab' && (!is.null(seurat_object.reactions$seurat@misc$cells_filtered) && seurat_object.reactions$seurat@misc$cells_filtered))
      sendSweetAlert(session=session, type='success', html=TRUE,
                     title='Notice', btn_labels='Great!',
                     text=tags$span('It looks like low-quality cells have already been removed from this Seurat object:', tags$h5(tags$code('@misc$cells_filtered == TRUE'))),
                     closeOnClickOutside=TRUE, showCloseButton=FALSE)
  })

  ## make knee plot of total expression
  cell_filtering.total_expression_knee.plot <- reactive(x={
    progress <- shiny::Progress$new(session=session, min=0, max=1/10)
    on.exit(progress$close())
    progress$set(value=0, message='Making total expression knee plot')

    min_value <- seurat_object.reactions$total_umi_per_cell_min
    max_value <- seurat_object.reactions$total_umi_per_cell_max

    progress$inc(detail='Making plot')
    FetchData(seurat_object.reactions$seurat, 'nCount_RNA') %>%
      set_names('y') %>%
      arrange(desc(y)) %>%
      mutate(x=seq(n()),
             pass=between(x=y, left=min_value, right=max_value)) %>%
      ggplot()+
      aes(x=x, y=y, colour=pass)+
      labs(x='Ranked cells', y='Total UMIs per cell', colour='Cells passing threshold')+
      geom_hline(yintercept=min_value, size=1) +
      geom_hline(yintercept=max_value, size=1) +
      geom_point(alpha=1, shape=16)+
      scale_x_log10(breaks=major_breaks_log10, minor_breaks=minor_breaks_log10, labels=function(x) scales::comma(x, accuracy=1))+
      scale_y_log10(breaks=major_breaks_log10, minor_breaks=minor_breaks_log10, labels=function(y) scales::comma(y, accuracy=1))+
      scale_colour_brewer(palette='Set1', direction=1) +
      theme_bw()+
      theme(legend.position='none', legend.title=element_blank())})

  ## make knee plot of unique genes detected
  cell_filtering.unique_genes_knee.plot <- reactive(x={
    progress <- shiny::Progress$new(session=session, min=0, max=1/10)
    on.exit(progress$close())
    progress$set(value=0, message='Making unique genes detected knee plot')

    min_value <- seurat_object.reactions$features_per_cell_min
    max_value <- seurat_object.reactions$features_per_cell_max

    progress$inc(detail='Making plot')
    FetchData(seurat_object.reactions$seurat, 'nFeature_RNA') %>%
      set_names('y') %>%
      arrange(desc(y)) %>%
      mutate(x=seq(n()),
             pass=between(x=y, left=min_value, right=max_value)) %>%
      ggplot()+
      aes(x=x, y=y, colour=pass)+
      labs(x='Ranked cells', y='Genes detected per cell', colour='Cells passing threshold')+
      geom_hline(yintercept=min_value, size=1) +
      geom_hline(yintercept=max_value, size=1) +
      geom_point(alpha=1, shape=16)+
      scale_x_log10(breaks=major_breaks_log10, minor_breaks=minor_breaks_log10, labels=function(x) scales::comma(x, accuracy=1))+
      scale_y_log10(breaks=major_breaks_log10, minor_breaks=minor_breaks_log10, labels=function(y) scales::comma(y, accuracy=1))+
      scale_colour_brewer(palette='Set1', direction=1) +
      theme_bw()+
      theme(legend.position='none', legend.title=element_blank())})

  ## make knee plot of mitochondrial expression
  cell_filtering.percent_mitochondria_knee.plot <- reactive(x={
    progress <- shiny::Progress$new(session=session, min=0, max=1/10)
    on.exit(progress$close())
    progress$set(value=0, message='Making mitochondrial expression knee plot')

    max_value <- seurat_object.reactions$percent_mt_per_cell_max

    progress$inc(detail='Making plot')
    # FetchData(seurat_object.reactions$seurat, 'percent_mt') %>%
    seurat_object.reactions$percent_mt %>%
      set_names('y') %>%
      arrange(desc(y)) %>%
      mutate(x=seq(n()),
             pass=y<=max_value) %>%
      ggplot()+
      aes(x=x, y=y, colour=pass)+
      labs(x='Ranked cells', y='Proportion mitochondrial expression', colour='Cells passing threshold')+
      geom_hline(yintercept=max_value, size=1) +
      geom_point(alpha=1, shape=16)+
      scale_x_log10(breaks=major_breaks_log10, minor_breaks=minor_breaks_log10, labels=function(x) scales::comma(x, accuracy=1))+
      scale_colour_brewer(palette='Set1', direction=1) +
      theme_bw()+
      theme(legend.position='none', legend.title=element_blank())})

  ## make density plot of total cell expression
  cell_filtering.total_expression_density.plot <- reactive(x={
    progress <- shiny::Progress$new(session=session, min=0, max=1/10)
    on.exit(progress$close())
    progress$set(value=0, message='Making total UMIs density plot')

    progress$inc(detail='Making plot')
    FetchData(seurat_object.reactions$seurat, 'nCount_RNA') %>%
      set_names('y') %>%
      ggplot() +
      aes(x=y) +
      labs(x='Total UMIs per cell', y='Density') +
      stat_density(geom='line', trim=TRUE, size=2) +
      scale_x_log10(breaks=major_breaks_log10, minor_breaks=minor_breaks_log10, labels=function(x) scales::comma(x, accuracy=1))+
      theme_bw()+
      theme()})

  ## make density plot of detected features
  cell_filtering.unique_features_density.plot <- reactive(x={
    progress <- shiny::Progress$new(session=session, min=0, max=1/10)
    on.exit(progress$close())
    progress$set(value=0, message='Making detected features density plot')

    progress$inc(detail='Making plot')
    FetchData(seurat_object.reactions$seurat, 'nFeature_RNA') %>%
      set_names('y') %>%
      ggplot() +
      aes(x=y) +
      labs(x='Detected features per cell', y='Density') +
      stat_density(geom='line', trim=TRUE, size=2) +
      scale_x_log10(breaks=major_breaks_log10, minor_breaks=minor_breaks_log10, labels=function(x) scales::comma(x, accuracy=1))+
      theme_bw()+
      theme()})

  ## make density plot of mitochondrial expression
  cell_filtering.percent_mitochondria_density.plot <- reactive(x={
    progress <- shiny::Progress$new(session=session, min=0, max=1/10)
    on.exit(progress$close())
    progress$set(value=0, message='Making percentage mitochondria density plot')

    progress$inc(detail='Making plot')
    # FetchData(seurat_object.reactions$seurat, 'percent_mt') %>%
    seurat_object.reactions$percent_mt %>%
      set_names('y') %>%
      ggplot() +
      aes(x=y) +
      labs(x='Proportion mitochondrial expression', y='Density') +
      stat_density(geom='line', trim=TRUE, size=2) +
      theme_bw()+
      theme()})

  ## make boxplot of total cell expression
  cell_filtering.total_expression_boxplot.plot <- reactive(x={
    progress <- shiny::Progress$new(session=session, min=0, max=1/10)
    on.exit(progress$close())
    progress$set(value=0, message='Making percentage mitochondria boxplot')

    min_y <- seurat_object.reactions$total_umi_per_cell_min
    max_y <- seurat_object.reactions$total_umi_per_cell_max

    progress$inc(detail='Making plot')
    FetchData(seurat_object.reactions$seurat, 'nCount_RNA') %>%
      set_names('y') %>%
      ggplot() +
      aes(x=0, y=y) +
      labs(y='Total UMIs per cell') +
      geom_point(shape=16, size=0.6, colour='lightgrey', alpha=0.6, position=position_jitter(width=0.5)) +
      geom_boxplot(fill=NA, width=0.4, size=0.9, outlier.size=0.6) +
      coord_cartesian(xlim=c(-1, 1), ylim=c(min_y, max_y)) +
      scale_y_continuous(labels=function(y) scales::comma(y, accuracy=1)) +
      theme_bw()+
      theme(axis.ticks.length.x=unit(0, 'pt'),
            axis.text.x=element_text(colour='white'),
            axis.title.x=element_text(colour='white'),
            panel.grid.major.x=element_blank(),
            panel.grid.minor.x=element_blank())})

  ## make boxplot of detected features
  cell_filtering.unique_genes_boxplot.plot <- reactive(x={
    progress <- shiny::Progress$new(session=session, min=0, max=1/10)
    on.exit(progress$close())
    progress$set(value=0, message='Making percentage mitochondria boxplot')

    min_y <- seurat_object.reactions$features_per_cell_min
    max_y <- seurat_object.reactions$features_per_cell_max

    progress$inc(detail='Making plot')
    FetchData(seurat_object.reactions$seurat, 'nFeature_RNA') %>%
      set_names('y') %>%
      ggplot() +
      aes(x=0, y=y) +
      labs(y='Detected features per cell') +
      geom_point(shape=16, size=0.6, colour='lightgrey', alpha=0.6, position=position_jitter(width=0.5)) +
      geom_boxplot(fill=NA, width=0.4, size=0.9, outlier.size=0.6) +
      coord_cartesian(xlim=c(-1, 1), ylim=c(min_y, max_y)) +
      scale_y_continuous(labels=function(y) scales::comma(y, accuracy=1)) +
      theme_bw()+
      theme(axis.ticks.length.x=unit(0, 'pt'),
            axis.text.x=element_text(colour='white'),
            axis.title.x=element_text(colour='white'),
            panel.grid.major.x=element_blank(),
            panel.grid.minor.x=element_blank())})

  ## make boxplot of mitochondrial expression
  cell_filtering.percent_mitochondria_boxplot.plot <- reactive(x={
    progress <- shiny::Progress$new(session=session, min=0, max=1/10)
    on.exit(progress$close())
    progress$set(value=0, message='Making percentage mitochondria boxplot')

    max_y <- seurat_object.reactions$percent_mt_per_cell_max

   progress$inc(detail='Making plot')
    # FetchData(seurat_object.reactions$seurat, 'percent_mt') %>%
    seurat_object.reactions$percent_mt %>%
      set_names('y') %>%
      ggplot() +
      aes(x=0, y=y) +
      labs(y='Proportion mitochondrial expression') +
      geom_point(shape=16, size=0.6, colour='lightgrey', alpha=0.6, position=position_jitter(width=0.5)) +
      geom_boxplot(fill=NA, width=0.4, size=0.9, outlier.size=0.6) +
      coord_cartesian(xlim=c(-1, 1), ylim=c(0, ceiling(max_y))) +
      theme_bw()+
      theme(axis.ticks.length.x=unit(0, 'pt'),
            axis.text.x=element_text(colour='white'),
            axis.title.x=element_text(colour='white'),
            panel.grid.major.x=element_blank(),
            panel.grid.minor.x=element_blank())})

  # ###############################################################################################
  # gene expression tab ---------------------------------------------------------------------------

  ## call modules
  ### when full palette type is selected
  callModule(module=update_palette_type.server, id='gene_highlighting')
  callModule(module=add_to_colour_palette.server, id='gene_highlighting')

  ## react to feature selection
  for(id in module_environments$feature_pickers$id)
    callModule(module=feature_picker.server, id=id)

  ## react to reduction method selection
  for(id in module_environments$reduction_method_pickers$id)
    callModule(module=reduction_method_picker.server, id=id)
  
  ## react to assay selection
  for(id in module_environments$assay_pickers$id)
    callModule(module=assay_picker.server, id=id)

  ## react to cluster set selection
  for(id in module_environments$cluster_resolution_pickers$id)
    callModule(module=cluster_resolution_picker.server, id=id)

  ## react to opacity selection
  for(id in module_environments$opacity_sliders$id)
    callModule(module=opacity_slider.server, id=id)

  ## react to point size selection
  for(id in module_environments$point_size_sliders$id)
    callModule(module=point_size_slider.server, id=id)

  ## gene highlighting
  # genes_highlighting.reactions <- reactiveValues(data=NULL, expression_map=NULL, cluster_expression=NULL, expression_map.running=0, cluster_expression.running=0)

  ## make map coloured by clusters
  genes_highlighting.clustered_map.plot <- reactive({
    progress <- shiny::Progress$new(session=session, min=0, max=if_else(seurat_object.reactions$label_clusters, 4, 2)/10)
    on.exit(progress$close())
    progress$set(value=0, message='Making cell clusters map')

    cluster_set <- sprintf('~%s', seurat_object.reactions$selected_cluster_resolution)

    progress$inc(detail='Making plot')
    cbind(seurat_object.reactions$dimred, seurat_object.reactions$picked_cluster_resolution_idents) %>%
      ggplot() +
      aes(x=DIMRED_1, y=DIMRED_2, colour=ident) +
      geom_hline(yintercept=0) + geom_vline(xintercept=0) +
      geom_point(size=seurat_object.reactions$point_size, alpha=seurat_object.reactions$opacity) +
      theme_void() +
      theme(legend.position='none') -> output_plot

    if(seurat_object.reactions$label_clusters) {
      progress$inc(detail='Getting cluster label positions')

      cbind(seurat_object.reactions$dimred, seurat_object.reactions$picked_cluster_resolution_idents) %>%
        group_by(ident) %>%
        summarise(DIMRED_1=mean(DIMRED_1), DIMRED_2=mean(DIMRED_2)) -> data_labels

      progress$inc(detail='Adding cluster labels')
      output_plot +
        ggrepel::geom_label_repel(data=data_labels,
                                  mapping=aes(label=ident),
                                  colour='black',
                                  size=12/(14/5)) -> output_plot
    }

    progress$inc(detail='Returning plot')
    output_plot})

  ## make map coloured by expression
  genes_highlighting.gene_expression_map.plot <- reactive({
    progress <- shiny::Progress$new(session=session, min=0, max=1/10)
    on.exit(progress$close())
    progress$set(value=0, message='Making expression per cluster plot')

    # update_slider()
    progress$inc(detail='Making plot')
    cbind(seurat_object.reactions$dimred, seurat_object.reactions$picked_feature_values) %>%
      rename(expression_value=value) -> data

    if(is.numeric(data$expression_value)) {
    data %>%
      arrange(expression_value) %>%
      ggplot() +
      aes(x=DIMRED_1, y=DIMRED_2, colour=expression_value) +
      geom_point(size=seurat_object.reactions$point_size, alpha=seurat_object.reactions$opacity) +
      scale_colour_gradient(low=input$`gene_highlighting-colour_palette-low`, high=input$`gene_highlighting-colour_palette-high`, limits=seurat_object.reactions$value_range_limits, oob=scales::squish) +
      # scale_colour_gradient(low=input$`gene_highlighting-colour_palette-low`, high=input$`gene_highlighting-colour_palette-high`, limits=input$expression_range.slider, oob=scales::squish) +
      theme_void() +
      theme(legend.position='none')
    } else {
    data %>%
      arrange(expression_value) %>%
      ggplot() +
      aes(x=DIMRED_1, y=DIMRED_2, colour=expression_value) +
      geom_point(size=seurat_object.reactions$point_size, alpha=seurat_object.reactions$opacity) +
      # scale_colour_gradient(low=input$expression_min.colour, high=input$expression_max.colour, limits=input$expression_range.slider, oob=scales::squish) +
      facet_wrap(~expression_value, scales='free') +
      guides(colour=guide_legend(override.aes=list(size=3, shape=15))) +
      theme_void() +
      theme(legend.position='bottom', legend.title=element_blank())
    }})

  ## plot expression ranges per cluster
  genes_highlighting.expression_per_cluster.plot <- reactive({
    progress <- shiny::Progress$new(session=session, min=0, max=3/10)
    on.exit(progress$close())
    progress$set(value=0, message='Making expression per cluster plot')

    # update_slider()
    progress$inc(detail='Fetching data')
    cbind(seurat_object.reactions$picked_cluster_resolution_idents, seurat_object.reactions$picked_feature_values) %>%
      rename(cluster_id=ident, expression_value=value) -> data

    if(is.numeric(data$expression_value)) {
      progress$inc(detail='Summarising expression in clusters')
      data %>%
        # filter(expression_value>0) %>%
        group_by(cluster_id) %>%
        summarise(q25=quantile(expression_value, 0.25), q75=quantile(expression_value, 0.75), median=median(expression_value)) %>%
        mutate(iqr=q75-q25, lower=q25-1.5*iqr, upper=q75+1.5*iqr) -> cluster_data_summary

      progress$inc(detail='Making plot')
      cluster_data_summary %>%
        gather(key='key', value='y', lower, upper) %>%
        mutate(x={as.character(cluster_id) %>% as.numeric()}) %>%
        ggplot() +
        aes(x=x, y=y, colour=cluster_id) +
        labs(x='Cluster identifier', y='Feature value (median ± 1.5x IQR)') +
        geom_line(size=1) +
        geom_point(mapping=aes(y=median), colour='black', shape=20, size=3) +
        theme_bw() +
        theme(legend.position='none', panel.grid.major.x=element_blank(), panel.grid.minor.x=element_blank())
    } else {
      data %>%
        mutate(x={as.character(cluster_id) %>% as.numeric()}) %>%
        ggplot() +
        aes(x=x, fill=expression_value) +
        labs(x='Cluster identifier', y='Frequency') +
        geom_bar(position='dodge') +
        theme_bw() +
        theme(legend.position='none', panel.grid.major.x=element_blank(), panel.grid.minor.x=element_blank())
    }})

  # ###############################################################################################
  # features heatmap tab --------------------------------------------------------------------------

  features_heatmap.heatmap.plot <- reactive(x={
    progress <- shiny::Progress$new(session=session, min=0, max=4/10)
    on.exit(progress$close())
    progress$set(value=0, message='Making heatmap of normalised feature values')

    # if(is.null(seurat_object.reactions$seurat))
      return(NULL)

    if(FALSE)
      input <- list(features_of_interest.tb={Matrix::rowSums(seurat_object.reactions$seurat) %>% (function(x) x[x>100]) %>% names()}, seurat_cluster_set.dd='integrated_snn_res.1.2')

    progress$inc(detail='Getting list of feature IDs')
    feature_names <- input$features_of_interest.tb
    if(feature_names[1]=='')
      Matrix::rowSums(seurat_object.reactions$seurat) %>%
        sort() %>%
        head(2000) %>%
        sample(size=200) %>%
        names() -> feature_names
    feature_names %<>% str_split(pattern='\\s|;|,') %>% unlist() %>% (function(x) x[x %in% rownames(seurat_object.reactions$seurat)])

    progress$inc(detail='Collecting feature values matrix')
    FetchData(object=seurat_object.reactions$seurat, vars=feature_names) %>%
      rownames_to_column('cell_id') -> feature_values_matrix
    FetchData(object=seurat_object.reactions$seurat, vars=input$features_heatmap.seurat_cluster_set.dd) %>%
      set_names('cluster_id') %>%
      rownames_to_column('cell_id') %>%
      arrange(cluster_id) %>%
      mutate(cell_n=seq(n())) -> cell_cluster_ids
    # cell_cluster_ids %>%
    #   group_by(cluster_id) %>%
    #   summarise(border_position=max(cell_n),
    #             label_position=mean(cell_n)) -> cluster_border_labels

    progress$inc(detail='Clustering features (columns)')
    # feature_values_matrix[,-1] %>%
    #   t() %>%
    #   dist() %>%
    #   hclust() %>%
    #   (function(x) {x$labels[x$order]}) -> ordered_feature_names

    # if(input$cluster_cells.checkbox) {
    # if(FALSE) {
    #   progress$inc(detail='Clustering cells (rows)')
    #   feature_values_matrix[,-1] %>%
    #     dist() %>%
    #     hclust() %>%
    #     (function(x) {x$labels[x$order]}) -> ordered_cell_ids
    # } else {
    #   ordered_cell_ids <- rownames(feature_values_matrix)
    # }

    progress$inc(detail='Making heatmap')
    # left_join(x=cell_cluster_ids, y=feature_values_matrix, by='cell_id') %>%
    #   gather(key='feature_name', value='value', -cell_id, -cell_n, -cluster_id) %>%
    #   mutate(feature_name=factor(x=feature_name, levels=rev(ordered_feature_names))) %>%
    #   ggplot() +
    #   aes(x=cell_n, y=feature_name, fill=value) +
    #   labs(x='Cluster ID', y='Feature name') +
    #   geom_raster() +
    #   geom_vline(data=cluster_border_labels, mapping=aes(xintercept=border_position), size=0.7, colour='grey') +
    #   scale_x_continuous(expand=c(0,0), labels=cluster_border_labels$cluster_id, breaks=cluster_border_labels$label_position, position='top') +
    #   scale_y_discrete(expand=c(0,0)) +
    #   scale_fill_viridis(limits=c(0,2), oob=scales::squish) +
    #   theme_bw() +
    #   theme(axis.ticks=element_blank(),
    #         legend.position='none',
    #         strip.background=element_blank()) -> heatmap_plot
    #
    # if(TRUE)
    #   heatmap_plot <- heatmap_plot + theme(axis.text.y=element_blank())

    if(!input$features_heatmap.use_average_value.checkbox) { # plot expression per gene per cell
      left_join(x=cell_cluster_ids, y=feature_values_matrix, by='cell_id') %>%
        select_at(vars('cluster_id', feature_names)) -> heatmap_data
    } else { # plot the mean expression per gene in cluster
    left_join(x=cell_cluster_ids, y=feature_values_matrix, by='cell_id') %>%
      group_by(cluster_id) %>%
      summarise_at(vars(feature_names), mean) -> heatmap_data
    }

    if(nrow(heatmap_data)<=1) # temporary fix
      return(NULL)

    heatmap_data %>%
      (function(x) Heatmap(matrix=t(scales::squish(x=as.matrix(x[,-1]), range=c(0,4))),
                           cluster_columns=FALSE,
                           cluster_rows=FALSE,
                           show_column_dend=FALSE,
                           show_row_dend=FALSE,
                           show_row_names=input$features_heatmap.show_feature_names.checkbox,
                           show_heatmap_legend=FALSE,
                           top_annotation=HeatmapAnnotation(`Cluster ID`=x$cluster_id, show_legend=FALSE),
                           col=viridis_pal()(30),
                           # height=unit(1, 'npc'),
                           width=unit(1, 'npc')))}) %>% debounce(2500)

  # ###############################################################################################
  # return plots to shiny -------------------------------------------------------------------------
  progress$inc(detail='Returning plots to Shiny')

  ## highlighting genes tab
  # renderPlotly({plot_ly(z=~volcano) %>% add_surface()}) -> output$`genes_highlighting-expression_per_cluster`

  callModule(module=project_name_text_box.server, id='gene_highlighting')
  callModule(module=number_of_cells_text_box.server, id='gene_highlighting')
  callModule(module=number_of_clusters_text_box.server, id='gene_highlighting')
  callModule(module=gene_name_and_description_text_box.server, id='gene_highlighting')
  callModule(module=number_of_reads_text_box.server, id='gene_highlighting')
  callModule(module=number_of_genes_in_assay_text_box.server, id='gene_highlighting')
  callModule(module=number_of_reads_per_cell_text_box.server, id='gene_highlighting')
  callModule(module=number_of_genes_per_cell_text_box.server, id='gene_highlighting')

  renderPlot(genes_highlighting.clustered_map.plot()) -> output$`genes_highlighting-clustered_map`
  renderPlot(genes_highlighting.gene_expression_map.plot()) -> output$`genes_highlighting-gene_expression_map`
  renderPlot(genes_highlighting.expression_per_cluster.plot()) -> output$`genes_highlighting-expression_per_cluster`

  ## cell filtering tab
  ### statistics boxes
  callModule(module=project_name_text_box.server, id='cell_filtering')
  callModule(module=number_of_reads_text_box.server, id='cell_filtering')
  callModule(module=number_of_cells_text_box.server, id='cell_filtering')
  callModule(module=number_of_reads_per_cell_text_box.server, id='cell_filtering')
  callModule(module=number_of_genes_per_cell_text_box.server, id='cell_filtering')

  ### knee plots
  renderPlot(cell_filtering.total_expression_knee.plot()) -> output$`cell_filtering-total_expression_knee`
  renderPlot(cell_filtering.unique_genes_knee.plot()) -> output$`cell_filtering-unique_genes_knee`
  renderPlot(cell_filtering.percent_mitochondria_knee.plot()) -> output$`cell_filtering-percent_mitochondria_knee`

  ### density plots
  renderPlot(cell_filtering.total_expression_density.plot()) -> output$`cell_filtering-total_expression_density`
  renderPlot(cell_filtering.unique_features_density.plot()) -> output$`cell_filtering-unique_genes_density`
  renderPlot(cell_filtering.percent_mitochondria_density.plot()) -> output$`cell_filtering-percent_mitochondria_density`

  ### boxplots
  renderPlot(cell_filtering.total_expression_boxplot.plot()) -> output$`cell_filtering-total_expression_boxplot`
  renderPlot(cell_filtering.unique_genes_boxplot.plot()) -> output$`cell_filtering-unique_genes_boxplot`
  renderPlot(cell_filtering.percent_mitochondria_boxplot.plot()) -> output$`cell_filtering-percent_mitochondria_boxplot`

  # features heatmap tab
  callModule(project_name_text_box.server, id='features_heatmap')
  renderPlot(features_heatmap.heatmap.plot()) -> output$`features_heatmap-heatmap`

  # sidebar
  ## dynamic sidebar outputs can be listed here

  # any code to exectue when the session ends
  session$onSessionEnded(function() {})
}
