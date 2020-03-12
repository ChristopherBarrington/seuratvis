
10^(0:9) -> major_breaks_log10
(2:9) * rep(major_breaks_log10, each=8) -> minor_breaks_log10

shinyAppServer <- function(input, output, session) {

  progress <- shiny::Progress$new(session=session, min=0, max=4/10)
  on.exit(progress$close())
  progress$set(value=0, message='Loading environment')

  # ###############################################################################################
  # load Seurat object from user ------------------------------------------------------------------

  ## react to Seurat object selection
  seurat_object.reactions <- reactiveValues()
  observeEvent(eventExpr=input$seurat_select.input, handlerExpr={

    progress <- shiny::Progress$new(session=session, min=0, max=7/10)
    on.exit(progress$close())
    progress$set(value=0, message='Loading environment')

    # load Seurat object from user
    progress$inc(detail='Updating UI with cluster options')
    seurat <- parse(text=input$seurat_select.input) %>% eval()
    seurat <- subset(seurat, subset=nFeature_RNA>0 & nCount_RNA>0)
    cluster_options <- c('seurat_clusters', str_subset(colnames(seurat@meta.data), '_snn_res.'))
    updateSelectInput(session=session, inputId='seurat_cluster_set.dd', choices=cluster_options)
    updateSelectInput(session=session, inputId='features_heatmap.seurat_cluster_set.dd', choices=cluster_options)

    progress$inc(detail='Combining UMAP and meta.data')
    if(!is.null(seurat@reductions$umap))
      umap <- cbind(seurat@meta.data, seurat@reductions$umap@cell.embeddings)

    if(is.null(seurat@meta.data$seurat_clusters))
      seurat@meta.data$seurat_clusters <- 0

    progress$inc(detail='Counting clusters identified in each set')
    select_at(seurat@meta.data, vars(contains('_snn_res.'), 'seurat_clusters')) %>%
      mutate_all(function(x) {as.character(x) %>% as.numeric()}) %>%
      gather(key='cluster_set', value='ID') %>%
      group_by(cluster_set) %>%
      summarise(N=length(unique(ID))+1) %>%
      deframe() -> clusters_per_resolution

    progress$inc(detail='Getting summary statistics')
    list(n_cells=nrow(seurat@meta.data),
         total_reads=sum(seurat@meta.data$nCount_RNA),
         median_reads_per_cell=round(x=median(seurat@meta.data$nCount_RNA), digits=0),
         median_genes_per_cell=round(x=median(seurat@meta.data$nFeature_RNA), digits=0),
         min_reads_per_cell=min(seurat@meta.data$nCount_RNA), max_reads_per_cell=max(seurat@meta.data$nCount_RNA),
         min_genes_per_cell=min(seurat@meta.data$nFeature_RNA), max_genes_per_cell=max(seurat@meta.data$nFeature_RNA),
         max_percent_mitochondria=round(max(seurat@meta.data$percent_mt)+0.05, digits=1)) -> cell_filtering_data.reference

    progress$inc(detail='Updating UI elements')
    updateTextInput(session=session, inputId='min_features_per_cell.textinput', placeholder=cell_filtering_data.reference$min_genes_per_cell)
    updateTextInput(session=session, inputId='max_features_per_cell.textinput', placeholder=cell_filtering_data.reference$max_genes_per_cell)
    updateTextInput(session=session, inputId='min_expression_per_cell.textinput', placeholder=cell_filtering_data.reference$min_reads_per_cell)
    updateTextInput(session=session, inputId='max_expression_per_cell.textinput', placeholder=cell_filtering_data.reference$max_reads_per_cell)
    updateTextInput(session=session, inputId='percent_mitochondria.textinput', placeholder=cell_filtering_data.reference$max_percent_mitochondria)
    update_autocomplete_input(session=session, id='gene_of_interest.dd', options=sort(rownames(seurat)))

    progress$inc(detail='Saving variables')
    available_assays <- Assays(seurat)
    available_slots <- lapply(seurat@assays, function(x) c('counts','data','scale.data') %>% purrr::set_names() %>% lapply(function(y) slot(x,y) %>% nrow())) %>% lapply(function(y) names(y)[unlist(y)>0])

    selected_assay <- 'RNA'
    selected_slot <- 'data'

    seurat@active.assay <- selected_assay

    # copy the important stuff into the reaction values
    seurat_object.reactions$seurat <- seurat
    seurat_object.reactions$mart <- seurat@misc$mart
    seurat_object.reactions$formatted.project.name <- seurat@project.name %>% str_replace_all(pattern='_', replacement=' ') %>% str_to_upper()
    seurat_object.reactions$reference_metrics <- cell_filtering_data.reference
    seurat_object.reactions$umap <- umap
    seurat_object.reactions$clusters_per_resolution <- clusters_per_resolution
  })

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
    min_genes_per_cell <- if_else(is.na(cell_filtering_data.reactions$min_genes_per_cell), as.numeric(seurat_object.reactions$reference_metrics$min_genes_per_cell), cell_filtering_data.reactions$min_genes_per_cell)
    max_genes_per_cell <- if_else(is.na(cell_filtering_data.reactions$max_genes_per_cell), as.numeric(seurat_object.reactions$reference_metrics$max_genes_per_cell), cell_filtering_data.reactions$max_genes_per_cell)
    min_expression_per_cell <- if_else(is.na(cell_filtering_data.reactions$min_expression_per_cell), as.numeric(seurat_object.reactions$reference_metrics$min_reads_per_cell), cell_filtering_data.reactions$min_expression_per_cell)
    max_expression_per_cell <- if_else(is.na(cell_filtering_data.reactions$max_expression_per_cell), as.numeric(seurat_object.reactions$reference_metrics$max_reads_per_cell), cell_filtering_data.reactions$max_expression_per_cell)
    max_percent_mitochondria <- if_else(is.na(cell_filtering_data.reactions$max_percent_mitochondria), as.numeric(seurat_object.reactions$reference_metrics$max_percent_mitochondria), cell_filtering_data.reactions$max_percent_mitochondria)

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

    # save formatted filters
    format_subset_conditional <- function(x, fmt) ifelse(is.na(x), NA, sprintf(fmt=fmt, x))
    group_format_subset_conditional <- function(x) x %>% na.omit() %>% paste(collapse=' & ')

    c(format_subset_conditional(x=min_expression_per_cell, fmt='nCount_RNA>=%d'),
      format_subset_conditional(x=max_expression_per_cell, fmt='nCount_RNA<=%d')) %>%
      group_format_subset_conditional() -> cell_filtering_data.reactions$subset_conditions$nCount

    c(format_subset_conditional(x=min_genes_per_cell, fmt='nFeature_RNA>=%d'),
      format_subset_conditional(x=max_genes_per_cell, fmt='nFeature_RNA<=%d')) %>%
      group_format_subset_conditional() -> cell_filtering_data.reactions$subset_conditions$nFeature

    format_subset_conditional(x=max_percent_mitochondria, fmt='percent_mt<=%s') -> cell_filtering_data.reactions$subset_conditions$percent_mt

    progress$inc(detail='Returning')
    NULL}) -> react_to_cell_filtering

  ## react to unique features detected density plot brush
  observeEvent(eventExpr=input$unique_genes_density.brush, handlerExpr={
    cell_filtering_data.reactions$min_genes_per_cell <- floor(input$unique_genes_density.brush$xmin)
    cell_filtering_data.reactions$max_genes_per_cell <- ceiling(input$unique_genes_density.brush$xmax)

    updateTextInput(session=session, inputId='min_features_per_cell.textinput', value=cell_filtering_data.reactions$min_genes_per_cell)
    updateTextInput(session=session, inputId='max_features_per_cell.textinput', value=cell_filtering_data.reactions$max_genes_per_cell)})

  observeEvent(eventExpr=input$min_features_per_cell.textinput, handlerExpr={
    # session$resetBrush(brushId='unique_genes_density.brush')
    cell_filtering_data.reactions$min_genes_per_cell <- floor(as.numeric(input$min_features_per_cell.textinput))})

  observeEvent(eventExpr=input$max_features_per_cell.textinput, handlerExpr={
    # session$resetBrush(brushId='unique_genes_density.brush')
    cell_filtering_data.reactions$max_genes_per_cell <- ceiling(as.numeric(input$max_features_per_cell.textinput))})

  ## react to total UMIs density plot brush
  observeEvent(eventExpr=input$total_expression_density.brush, handlerExpr={
    cell_filtering_data.reactions$min_expression_per_cell <- floor(input$total_expression_density.brush$xmin)
    cell_filtering_data.reactions$max_expression_per_cell <- ceiling(input$total_expression_density.brush$xmax)

    updateTextInput(session=session, inputId='min_expression_per_cell.textinput', value=cell_filtering_data.reactions$min_expression_per_cell)
    updateTextInput(session=session, inputId='max_expression_per_cell.textinput', value=cell_filtering_data.reactions$max_expression_per_cell)})

  observeEvent(eventExpr=input$min_expression_per_cell.textinput, handlerExpr={
    # session$resetBrush(brushId='total_expression_density.brush')
    cell_filtering_data.reactions$min_expression_per_cell <- floor(as.numeric(input$min_expression_per_cell.textinput))})

  observeEvent(eventExpr=input$max_expression_per_cell.textinput, handlerExpr={
    # session$resetBrush(brushId='total_expression_density.brush')
    cell_filtering_data.reactions$max_expression_per_cell <- ceiling(as.numeric(input$max_expression_per_cell.textinput))})

  ## react to percent mitochondria density plot brush
  observeEvent(eventExpr=input$percent_mitochondria_density.brush, handlerExpr={
    cell_filtering_data.reactions$max_percent_mitochondria <- round(input$percent_mitochondria_density.brush$xmax+0.05, digits=1)

    updateTextInput(session=session, inputId='percent_mitochondria.textinput', value=cell_filtering_data.reactions$max_percent_mitochondria)})

  observeEvent(eventExpr=input$percent_mitochondria.textinput, handlerExpr={
    # session$resetBrush(brushId='percent_mitochondria_density.brush')
    cell_filtering_data.reactions$max_percent_mitochondria <- as.numeric(input$percent_mitochondria.textinput)})

  observeEvent(input$sidebarmenu, {
    if(input$sidebarmenu=='cell_filtering-tab' & (!is.null(seurat_object.reactions$seurat@misc$cells_filtered) && seurat_object.reactions$seurat@misc$cells_filtered))
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

    min_value <- if_else(is.na(cell_filtering_data.reactions$min_expression_per_cell), seurat_object.reactions$reference_metrics$min_reads_per_cell, cell_filtering_data.reactions$min_expression_per_cell)
    max_value <- if_else(is.na(cell_filtering_data.reactions$max_expression_per_cell), seurat_object.reactions$reference_metrics$max_reads_per_cell, cell_filtering_data.reactions$max_expression_per_cell)

    progress$inc(detail='Making plot')
    FetchData(seurat_object.reactions$seurat, 'nCount_RNA') %>%
      set_names('y') %>%
      arrange(desc(y)) %>%
      mutate(x=seq(n()),
             pass=between(x=y, left=min_value, right=max_value)) %>%
      ggplot()+
      aes(x=x, y=y, colour=pass)+
      labs(x='Ranked cells', y='Total expression per cell', colour='Cells passing threshold')+
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

    min_value <- if_else(is.na(cell_filtering_data.reactions$min_genes_per_cell), as.numeric(seurat_object.reactions$reference_metrics$min_genes_per_cell), cell_filtering_data.reactions$min_genes_per_cell)
    max_value <- if_else(is.na(cell_filtering_data.reactions$max_genes_per_cell), as.numeric(seurat_object.reactions$reference_metrics$max_genes_per_cell), cell_filtering_data.reactions$max_genes_per_cell)

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

    max_value <- if_else(is.na(cell_filtering_data.reactions$max_percent_mitochondria), seurat_object.reactions$reference_metrics$max_percent_mitochondria, cell_filtering_data.reactions$max_percent_mitochondria)

    progress$inc(detail='Making plot')
    FetchData(seurat_object.reactions$seurat, 'percent_mt') %>%
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
    FetchData(seurat_object.reactions$seurat, 'percent_mt') %>%
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

    min_y <- if_else(is.na(cell_filtering_data.reactions$min_expression_per_cell), seurat_object.reactions$reference_metrics$min_reads_per_cell, cell_filtering_data.reactions$min_expression_per_cell)
    max_y <- if_else(is.na(cell_filtering_data.reactions$max_expression_per_cell), seurat_object.reactions$reference_metrics$max_reads_per_cell, cell_filtering_data.reactions$max_expression_per_cell)

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

    min_y <- if_else(is.na(cell_filtering_data.reactions$min_genes_per_cell), as.numeric(seurat_object.reactions$reference_metrics$min_genes_per_cell), cell_filtering_data.reactions$min_genes_per_cell)
    max_y <- if_else(is.na(cell_filtering_data.reactions$max_genes_per_cell), as.numeric(seurat_object.reactions$reference_metrics$max_genes_per_cell), cell_filtering_data.reactions$max_genes_per_cell)

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

    max_y <- if_else(is.na(cell_filtering_data.reactions$max_percent_mitochondria), seurat_object.reactions$reference_metrics$max_percent_mitochondria, cell_filtering_data.reactions$max_percent_mitochondria)

   progress$inc(detail='Making plot')
    FetchData(seurat_object.reactions$seurat, 'percent_mt') %>%
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

  ## update UI
  ### when gene is selected
  reactive(x={
    progress <- shiny::Progress$new(session=session, min=0, max=1/10)
    on.exit(progress$close())
    progress$set(value=0, message=sprintf('Updating UI for %s', input$gene_of_interest.dd))

    # cat(file=stderr(), sprintf('# reacting to gene selection with: %s\n', input$gene_of_interest.dd))
    updateSliderInput(session=session, inputId='expression_range.slider', max={max(FetchData(object=seurat_object.reactions$seurat, vars=input$gene_of_interest.dd)+0.05) %>% round(digits=1)}, value=c(-Inf,Inf))

    progress$inc(detail='Done')
    NULL}) -> update_slider

  ### when palette full palette type is selected
  observeEvent(eventExpr=input$expression_palette_full, handlerExpr={
    progress <- shiny::Progress$new(session=session, min=0, max=1/10)
    on.exit(progress$close())
    progress$set(value=0, message='Updating colour palette UI')

    palette_type <- ifelse(input$expression_palette_full, 'square', 'limited')

    for(inputId in c('expression_min.colour','expression_max.colour'))
      updateColourInput(session=session, inputId=inputId, palette=palette_type, value=input[[inputId]], allowedCols=colour_palette)

    progress$inc(detail='Done')}) -> update_palette_type

  ### add new colours to palette
  observeEvent(eventExpr=input$expression_min.colour, handlerExpr={
    add_to_colour_palette(input$expression_min.colour)})

  observeEvent(eventExpr=input$expression_max.colour, handlerExpr={
    add_to_colour_palette(input$expression_max.colour)})

  ## gene highlighting
  # genes_highlighting.reactions <- reactiveValues(data=NULL, expression_map=NULL, cluster_expression=NULL, expression_map.running=0, cluster_expression.running=0)

  ## make map coloured by clusters
  genes_highlighting.clustered_map.plot <- reactive({
    progress <- shiny::Progress$new(session=session, min=0, max=if_else(input$gene_highlighting.label_clusters.checkbox, 4, 2)/10)
    on.exit(progress$close())
    progress$set(value=0, message='Making cell clusters map')

    cluster_set <- sprintf('~%s', input$seurat_cluster_set.dd)

    progress$inc(detail='Making plot')
    seurat_object.reactions$umap %>%
      ggplot() +
      aes(x=UMAP_1, y=UMAP_2) +
      aes_string(colour=input$seurat_cluster_set.dd) +
      geom_point(size=input$gene_highlighting.point_size.slider, alpha=input$opacity.slider) +
      theme_void() +
      theme(legend.position='none') -> output_plot

    if(input$gene_highlighting.label_clusters.checkbox) {
      progress$inc(detail='Getting cluster label positions')
      seurat_object.reactions$umap %>%
        group_by_at(vars(cluster_id=input$seurat_cluster_set.dd)) %>%
        summarise(UMAP_1=mean(UMAP_1), UMAP_2=mean(UMAP_2)) -> data_labels

      progress$inc(detail='Adding cluster labels')
      output_plot +
        ggrepel::geom_label_repel(data=data_labels,
                                  mapping=aes(label=cluster_id),
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

    update_slider()
    progress$inc(detail='Making plot')
    FetchData(object=seurat_object.reactions$seurat, vars=c('UMAP_1','UMAP_2',input$gene_of_interest.dd)) %>%
      set_names(c('UMAP_1','UMAP_2', 'expression_value')) %>%
      arrange(expression_value) %>%
      ggplot() +
      aes(x=UMAP_1, y=UMAP_2, colour=expression_value) +
      geom_point(size=input$gene_highlighting.point_size.slider, alpha=input$opacity.slider) +
      scale_colour_gradient(low=input$expression_min.colour, high=input$expression_max.colour, limits=input$expression_range.slider, oob=scales::squish) +
      theme_void() +
      theme(legend.position='none')})

  ## plot expression ranges per cluster
  genes_highlighting.expression_per_cluster.plot <- reactive({
    progress <- shiny::Progress$new(session=session, min=0, max=3/10)
    on.exit(progress$close())
    progress$set(value=0, message='Making expression per cluster plot')

    update_slider()
    progress$inc(detail='Fetching data')
    FetchData(object=seurat_object.reactions$seurat, vars=c(input$seurat_cluster_set.dd,input$gene_of_interest.dd)) %>%
      set_names(c('cluster_id', 'expression_value')) -> data

    progress$inc(detail='Summarising expression in clusters')
    data %>%
      group_by(cluster_id) %>%
      summarise(q25=quantile(expression_value, 0.25), q75=quantile(expression_value, 0.75), median=median(expression_value)) %>%
      mutate(iqr=q75-q25, lower=q25-1.5*iqr, upper=q75+1.5*iqr) -> cluster_data_summary

    progress$inc(detail='Making plot')
    cluster_data_summary %>%
      gather(key='key', value='y', lower, upper) %>%
      mutate(x={as.character(cluster_id) %>% as.numeric()}) %>%
      ggplot() +
      aes(x=x, y=y, colour=cluster_id) +
      labs(x='Cluster identifier', y='Normalised expression') +
      geom_line(size=1) +
      geom_point(mapping=aes(y=median), colour='black', shape=20, size=3) +
      theme_bw() +
      theme(legend.position='none', panel.grid.major.x=element_blank(), panel.grid.minor.x=element_blank())})

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
  renderPlot(genes_highlighting.clustered_map.plot()) -> output$`genes_highlighting-clustered_map`
  renderPlot(genes_highlighting.gene_expression_map.plot()) -> output$`genes_highlighting-gene_expression_map`
  renderPlot(genes_highlighting.expression_per_cluster.plot()) -> output$`genes_highlighting-expression_per_cluster`

  output$genes_highlighting.n_cells_box <- renderValueBox(expr={
    valueBox(value=scales::comma(ncol(seurat_object.reactions$seurat)),
             subtitle='Cells in map',
             icon=icon('galactic-republic'),
             color='purple')})
  output$genes_highlighting.n_clusters_box <- renderValueBox(expr={
    valueBox(value=scales::comma(seurat_object.reactions$clusters_per_resolution[input$seurat_cluster_set.dd]),
             subtitle='Cell clusters',
             icon=icon('first-order'),
             color='purple')})
  output$genes_highlighting.selected_gene_box <- renderValueBox(expr={
    valueBox(value=input$gene_of_interest.dd,
             subtitle=get_formatted_gene_description(mart=seurat_object.reactions$mart, external_gene_name=input$gene_of_interest.dd),
             icon=icon('jedi-order'),
             color='purple')})
  output$genes_highlighting.n_genes_box <- renderValueBox(expr={
    valueBox(value=scales::comma(sapply(seurat_object.reactions$seurat@assays, function(x) nrow(x@data)) %>% max()),
             subtitle='Unique genes detected (max)',
             icon=icon('galactic-senate'),
             color='purple')})
  output$genes_highlighting.n_reads_box <- renderValueBox(expr={
    valueBox(value=scales::comma(sum(seurat_object.reactions$seurat$nCount_RNA)),
             subtitle='Total reads in cells',
             icon=icon('old-republic'),
             color='purple')})
  output$genes_highlighting.n_reads_per_cell_box <- renderValueBox(expr={
    valueBox(value=scales::comma(round(median(seurat_object.reactions$seurat$nCount_RNA), digits=1)),
             subtitle='Median reads per cell',
             icon=icon('frog'),
             color='purple')})
  output$genes_highlighting.n_genes_per_cell_box <- renderValueBox(expr={
    valueBox(value=scales::comma(round(median(seurat_object.reactions$seurat$nFeature_RNA), digits=1)),
             subtitle='Median genes per cell',
             icon=icon('crow'),
             color='purple')})

  ## cell filtering tab
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

  ### statistics boxes
  output$cell_filtering.n_reads_box <- renderValueBox({
    react_to_cell_filtering()
    valueBox(value=scales::comma(cell_filtering_data.reactions$total_reads),
             subtitle=sprintf(fmt='Total reads remaining (%.1f%%)', cell_filtering_data.reactions$total_reads/seurat_object.reactions$reference_metrics$total_reads*100),
             icon=icon('old-republic'),
             color='purple')})
  output$cell_filtering.n_cells_box <- renderValueBox({
    react_to_cell_filtering()
    valueBox(value=scales::comma(cell_filtering_data.reactions$n_cells),
             subtitle=sprintf(fmt='Cells remaining (%.1f%%)', cell_filtering_data.reactions$n_cells/seurat_object.reactions$reference_metrics$n_cells*100),
             icon=icon('galactic-republic'),
             color='purple')})
  output$cell_filtering.n_reads_per_cell_box <- renderValueBox({
    react_to_cell_filtering()
    valueBox(value=scales::comma(cell_filtering_data.reactions$median_reads_per_cell),
             subtitle=sprintf(fmt='Median reads per cell (%s)', comma(cell_filtering_data.reactions$median_reads_per_cell-seurat_object.reactions$reference_metrics$median_reads_per_cell) %>% ifelse(str_detect(., '^-'), ., str_c('+', .))),
             icon=icon('frog'),
             color='purple')})
  output$cell_filtering.n_genes_per_cell_box <- renderValueBox({
    react_to_cell_filtering()
    valueBox(value=scales::comma(cell_filtering_data.reactions$median_genes_per_cell),
             subtitle=sprintf(fmt='Median genes per cell (%s)', comma(cell_filtering_data.reactions$median_genes_per_cell-seurat_object.reactions$reference_metrics$median_genes_per_cell) %>% ifelse(str_detect(., '^-'), ., str_c('+', .))),
             icon=icon('crow'),
             color='purple')})

  ### formatted text box with filtering parameters
  output$`cell_filtering-subset_conditions` <- renderText({
    react_to_cell_filtering()
    c(sprintf(fmt='# %s', seurat_object.reactions$seurat@project.name),
      sprintf(fmt='# n_cells=%s', comma(cell_filtering_data.reactions$n_cells)),
      paste(cell_filtering_data.reactions$subset_conditions, collapse=' &\n')) %>%
      paste(collapse='\n') -> cell_filtering_data.reactions$subset_conditions.text})

  ### copy to clipboard buttons
  output$`cell_filtering-subset_conditions.plain` <- renderUI({ # copy text as is
    react_to_cell_filtering()
    rclipButton(inputId='rclipButton.plain.in', label='', icon('clipboard-check'),
                clipText=cell_filtering_data.reactions$subset_conditions.text)})

  output$`cell_filtering-subset_conditions.tsv` <- renderUI({ # tab separate elements
    react_to_cell_filtering()
    c(seurat_object.reactions$seurat@project.name,
      cell_filtering_data.reactions$n_cells,
      paste(cell_filtering_data.reactions$subset_conditions, collapse=' & ')) %>%
      paste(collapse=',') %>%
      rclipButton(inputId='rclipButton.tsv.in', label='', icon('file-excel'))})

  output$`cell_filtering-subset_conditions.r` <- renderUI({ # copy only the conditional
    react_to_cell_filtering()
    paste(cell_filtering_data.reactions$subset_conditions, collapse=' & ') %>%
      rclipButton(inputId='rclipButton.r.in', label='', icon('r-project'))})

  # features heatmap tab
  renderPlot(features_heatmap.heatmap.plot()) -> output$`features_heatmap-heatmap`

  # sidebar
  ## dynamic sidebar outputs can be listed here

  # shared text boxes

  ## project name
  project_name_box_opts <- list(subtitle='Loaded Seurat object', icon=icon('certificate'), color='purple')
  output$cell_filtering.project_name_box <- renderValueBox(expr={append(project_name_box_opts, list(value=seurat_object.reactions$formatted.project.name)) %>% do.call(what=valueBox)})
  output$genes_highlighting.project_name_box <- renderValueBox(expr={append(project_name_box_opts, list(value=seurat_object.reactions$formatted.project.name)) %>% do.call(what=valueBox)})
  output$features_heatmap.project_name_box <- renderValueBox(expr={append(project_name_box_opts, list(value=seurat_object.reactions$formatted.project.name)) %>% do.call(what=valueBox)})

  # any code to exectue when the session ends
  session$onSessionEnded(function() {
    colour_palette <<- default_colour_palette() # do not carry custom palettes between sessions
  })
}
