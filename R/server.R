
# library(promises)
# library(future)
# plan(multiprocess)

# library(biomaRt) ## move this into the seurat object@results slot

# library(hdf5r)
# library(shiny)
# library(ggvis)
# library(viridis)
# library(ComplexHeatmap)
# library(tidyverse)
# library(magrittr)

# load('int.RData') ######
# # seurat <- human_CS17_thoracic
# seurat <- human_CS17_brachial

read.delim('D:/seurat-vis/temp-gene-info.tsv') %>%
  filter(description!='') %>%
  mutate(description=str_remove_all(description, pattern=' \\[Source.+$')) %>%
  dplyr::select(external_gene_name, description) %>%
  unique() %>%
  deframe() -> gene_names_to_description

10^(0:9) -> major_breaks_log10
(2:9) * rep(major_breaks_log10, each=8) -> minor_breaks_log10

shinyAppServer <- function(input, output, session) {

#####################################################################################################################
# this is just placeholder stuff
  output$plot1 <- renderPlot({

    plot1.data <<- rnorm(input$slider)
    hist(plot1.data)
  })

  output$progressBox <- renderValueBox({
    valueBox(value=sprintf('%.1f', sum(plot1.data>0)/input$slider*100),
             subtitle='Values above 0',
             icon=icon('percent'),
             color="purple")})
  
  output$plot2 <- renderPlot({
    data <- runif(input$slider2)
    hist(data)})
  
  output$menuitem <- renderMenu({
    s <- stringi::stri_rand_strings(n=1, length=5)
    s <- Sys.getenv('USER')
    menuItem(text=s, icon=icon("file-code-o"), href="https://github.com/rstudio/shinydashboard")})
# this was just placeholder stuff
#####################################################################################################################

  progress <- shiny::Progress$new(session=session, min=0, max=4/10)
  on.exit(progress$close())
  progress$set(value=0, message='Loading environment')

  # load Seurat object from user
  progress$inc(detail='Updating UI with cluster options')
  cluster_options <- c('seurat_clusters', str_subset(colnames(seurat@meta.data), '_snn_res.'))
  updateSelectInput(session=session, inputId='seurat_cluster_set.dd', choices=cluster_options)
  updateSelectInput(session=session, inputId='features_heatmap.seurat_cluster_set.dd', choices=cluster_options)
  
  progress$inc(detail='Combining UMAP and meta.data')
  if(!is.null(seurat@reductions$umap))
    seurat_metadata_umap <- cbind(seurat@meta.data, seurat@reductions$umap@cell.embeddings)

  progress$inc(detail='Counting clusters identified in each set')
  select_at(seurat@meta.data, vars(contains('_snn_res.'), 'seurat_clusters')) %>%
    mutate_all(function(x) {as.character(x) %>% as.numeric()}) %>%
    gather(key='cluster_set', value='ID') %>%
    group_by(cluster_set) %>%
    summarise(N=length(unique(ID))+1) %>%
    deframe() -> seurat_cluster_per_set

  list(n_cells=nrow(seurat@meta.data),
       total_reads=sum(seurat@meta.data$nCount_RNA),
       median_reads_per_cell=round(x=median(seurat@meta.data$nCount_RNA), digits=0),
       median_genes_per_cell=round(x=median(seurat@meta.data$nFeature_RNA), digits=0)) -> cell_filtering_data.reference

  available_assays <- Assays(seurat)
  available_slots <- lapply(seurat@assays, function(x) c('counts','data','scale.data') %>% purrr::set_names() %>% lapply(function(y) slot(x,y) %>% nrow())) %>% lapply(function(y) names(y)[unlist(y)>0])
  
  selected_assay <- 'RNA'
  selected_slot <- 'data'
  
  seurat@active.assay <- selected_assay
  
  # ###############################################################################################
  # cell filtering tab ----------------------------------------------------------------------------

  ## react to cell filtering parameters
  cell_filtering_data.reactions <- reactiveValues(filtered_cell_set=NULL, n_cells=0, total_reads=0, median_reads_per_cell=0, median_genes_per_cell=0)
  reactive(x={
    progress <- shiny::Progress$new(session=session, min=0, max=4/10)
    on.exit(progress$close())
    progress$set(value=0, message='Reacting to threshold parameters')

    progress$inc(detail='Extracting input options')
    min_genes_per_cell <- input$min_genes_per_cell.slider
    max_genes_per_cell <- input$max_genes_per_cell.slider
    min_expression_per_cell <- input$min_expression_per_cell.slider
    max_expression_per_cell <- input$max_expression_per_cell.slider

    progress$inc(detail='Filtering @meta.data')
    seurat@meta.data %>%
      filter(percent_mt <= input$percent_mitochondria.slider &
             between(x=nFeature_RNA, left=min_genes_per_cell, right=max_genes_per_cell) &
             between(x=nCount_RNA, left=min_expression_per_cell, right=max_expression_per_cell)) -> filtered_cell_set

    progress$inc(detail='Saving results to reactiveValues')
    cell_filtering_data.reactions$filtered_cell_set <- filtered_cell_set
    cell_filtering_data.reactions$n_cells <- nrow(filtered_cell_set)
    cell_filtering_data.reactions$total_reads <- sum(filtered_cell_set$nCount_RNA)
    cell_filtering_data.reactions$median_reads_per_cell <- round(x=median(filtered_cell_set$nCount_RNA), digits=0)
    cell_filtering_data.reactions$median_genes_per_cell <- round(x=median(filtered_cell_set$nFeature_RNA), digits=0)

    progress$inc(detail='Returning')
    NULL}) -> react_to_cell_filtering

  ## make knee plot of total expression
  cell_filtering.total_expression_knee.plot <- reactive(x={
    progress <- shiny::Progress$new(session=session, min=0, max=1/10)
    on.exit(progress$close())
    progress$set(value=0, message='Making total expression knee plot')

    progress$inc(detail='Making plot')
    FetchData(seurat, 'nCount_RNA') %>%
      set_names('y') %>%
      arrange(desc(y)) %>%
      mutate(x=seq(n()),
             pass=between(x=y, left=input$min_expression_per_cell.slider, right=input$max_expression_per_cell.slider)) %>%
      ggplot()+
      aes(x=x, y=y, colour=pass)+
      labs(x='Ranked cells', y='Total expression per cell', colour='Cells passing threshold')+
      geom_hline(yintercept=input$min_expression_per_cell.slider, size=1) +
      geom_hline(yintercept=input$max_expression_per_cell.slider, size=1) +
      geom_point(alpha=1, shape=16)+
      scale_y_log10(minor_breaks=minor_breaks_log10, labels=scales::comma)+
      scale_x_log10(minor_breaks=minor_breaks_log10, labels=scales::comma)+
      scale_colour_brewer(palette='Set1', direction=1) +
      theme_bw()+
      theme(legend.position='none', legend.title=element_blank())})

  ## make knee plot of unique genes detected
  cell_filtering.unique_genes_knee.plot <- reactive(x={
    progress <- shiny::Progress$new(session=session, min=0, max=1/10)
    on.exit(progress$close())
    progress$set(value=0, message='Making unique genes detected knee plot')

    progress$inc(detail='Making plot')
    FetchData(seurat, 'nFeature_RNA') %>%
      set_names('y') %>%
      arrange(desc(y)) %>%
      mutate(x=seq(n()),
             pass=between(x=y, left=input$min_genes_per_cell.slider, right=input$max_genes_per_cell.slider)) %>%
      ggplot()+
      aes(x=x, y=y, colour=pass)+
      labs(x='Ranked cells', y='Genes detected per cell', colour='Cells passing threshold')+
      geom_hline(yintercept=input$min_genes_per_cell.slider, size=1) +
      geom_hline(yintercept=input$max_genes_per_cell.slider, size=1) +
      geom_point(alpha=1, shape=16)+
      scale_y_log10(minor_breaks=minor_breaks_log10, labels=scales::comma)+
      scale_x_log10(minor_breaks=minor_breaks_log10, labels=scales::comma)+
      scale_colour_brewer(palette='Set1', direction=1) +
      theme_bw()+
      theme(legend.position='none', legend.title=element_blank())})

  ## make knee plot of mitochondrial expression
  cell_filtering.percent_mitochondria_knee.plot <- reactive(x={
    progress <- shiny::Progress$new(session=session, min=0, max=1/10)
    on.exit(progress$close())
    progress$set(value=0, message='Making mitochondrial expression knee plot')
  
    progress$inc(detail='Making plot')
    FetchData(seurat, 'percent_mt') %>%
      set_names('y') %>%
      arrange(desc(y)) %>%
      mutate(x=seq(n()),
             pass=y<=input$percent_mitochondria.slider) %>%
      ggplot()+
      aes(x=x, y=y, colour=pass)+
      labs(x='Ranked cells', y='Proportion mitochondrial expression', colour='Cells passing threshold')+
      geom_hline(yintercept=input$percent_mitochondria.slider, size=1) +
      geom_point(alpha=1, shape=16)+
      scale_x_log10(minor_breaks=minor_breaks_log10, labels=scales::comma)+
      scale_colour_brewer(palette='Set1', direction=1) +
      theme_bw()+
      theme(legend.position='none', legend.title=element_blank())})

  ## make density plot of total cell expression
  cell_filtering.total_expression_density.plot <- reactive(x={
    progress <- shiny::Progress$new(session=session, min=0, max=1/10)
    on.exit(progress$close())
    progress$set(value=0, message='Making total UMIs density plot')
  
    progress$inc(detail='Making plot')
    FetchData(seurat, 'nCount_RNA') %>%
      set_names('y') %>%
      ggplot() +
      aes(x=y) +
      labs(x='Total UMIs per cell') +
      geom_vline(xintercept=input$min_expression_per_cell.slider, colour='#377eb8', size=1) +
      geom_vline(xintercept=input$max_expression_per_cell.slider, colour='#e41a1c', size=1) +
      stat_density(geom='line', trim=TRUE, size=2) +
      scale_x_log10(minor_breaks=minor_breaks_log10, labels=scales::comma)+
      theme_bw()+
      theme()})

  ## make density plot of detected features
  cell_filtering.unique_features_density.plot <- reactive(x={
    progress <- shiny::Progress$new(session=session, min=0, max=1/10)
    on.exit(progress$close())
    progress$set(value=0, message='Making detected features density plot')
  
    progress$inc(detail='Making plot')
    FetchData(seurat, 'nFeature_RNA') %>%
      set_names('y') %>%
      ggplot() +
      aes(x=y) +
      labs(x='Detected features per cell') +
      geom_vline(xintercept=input$min_genes_per_cell.slider, colour='#377eb8', size=1) +
      geom_vline(xintercept=input$max_genes_per_cell.slider, colour='#e41a1c', size=1) +
      stat_density(geom='line', trim=TRUE, size=2) +
      scale_x_log10(minor_breaks=minor_breaks_log10, labels=scales::comma)+
      theme_bw()+
      theme()})

  ## make density plot of detected features
  cell_filtering.percent_mitochondria_density.plot <- reactive(x={
    progress <- shiny::Progress$new(session=session, min=0, max=1/10)
    on.exit(progress$close())
    progress$set(value=0, message='Making percentage mitochondria density plot')
  
    progress$inc(detail='Making plot')
    FetchData(seurat, 'percent_mt') %>%
      set_names('y') %>%
      ggplot() +
      aes(x=y) +
      labs(x='Proportion mitochondrial expression') +
      geom_vline(xintercept=input$percent_mitochondria.slider, colour='#e41a1c', size=1) +
      stat_density(geom='line', trim=TRUE, size=2) +
      theme_bw()+
      theme()})

  # ###############################################################################################
  # gene expression tab ---------------------------------------------------------------------------
  
  ## update UI when gene is selected
  reactive(x={
    progress <- shiny::Progress$new(session=session, min=0, max=1/10)
    on.exit(progress$close())
    progress$set(value=0, message=sprintf('Updating UI for %s', input$gene_of_interest.dd))

    # cat(file=stderr(), sprintf('# reacting to gene selection with: %s\n', input$gene_of_interest.dd))
    updateSliderInput(session=session, inputId='expression_range.slider', max={max(FetchData(object=seurat, vars=input$gene_of_interest.dd)+0.05) %>% round(digits=1)}, value=c(-Inf,Inf))
    
    progress$inc(detail='Done')
    NULL}) -> update_slider

  ## gene highlighting
  # genes_highlighting.reactions <- reactiveValues(data=NULL, expression_map=NULL, cluster_expression=NULL, expression_map.running=0, cluster_expression.running=0)

  ## make map coloured by clusters
  genes_highlighting.clustered_map.plot <- reactive({
    progress <- shiny::Progress$new(session=session, min=0, max=if_else(input$gene_highlighting.label_clusters.checkbox, 4, 2)/10)
    on.exit(progress$close())
    progress$set(value=0, message='Making cell clusters map')

    cluster_set <- sprintf('~%s', input$seurat_cluster_set.dd)

    progress$inc(detail='Making plot')
    seurat_metadata_umap %>%
      ggplot() +
      aes(x=UMAP_1, y=UMAP_2) +
      aes_string(colour=input$seurat_cluster_set.dd) +
      geom_point(size=input$gene_highlighting.point_size.slider, alpha=input$opacity.slider) +
      theme_void() +
      theme(legend.position='none') -> output_plot

    if(input$gene_highlighting.label_clusters.checkbox) {
      progress$inc(detail='Getting cluster label positions')
      seurat_metadata_umap %>%
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
    FetchData(object=seurat, vars=c('UMAP_1','UMAP_2',input$gene_of_interest.dd)) %>%
      set_names(c('UMAP_1','UMAP_2', 'expression_value')) %>%
      arrange(expression_value) %>%
      ggplot() +
      aes(x=UMAP_1, y=UMAP_2, colour=expression_value) +
      geom_point(size=input$gene_highlighting.point_size.slider, alpha=input$opacity.slider) +
      scale_colour_gradient(low=input$expression_min.colour.dd, high=input$expression_max.colour.dd, limits=input$expression_range.slider, oob=scales::squish) +
      theme_void() +
      theme(legend.position='none')})
  
  ## plot expression ranges per cluster
  genes_highlighting.expression_per_cluster.plot <- reactive({
    progress <- shiny::Progress$new(session=session, min=0, max=3/10)
    on.exit(progress$close())
    progress$set(value=0, message='Making expression per cluster plot')

    update_slider()
    progress$inc(detail='Fetching data')
    FetchData(object=seurat, vars=c(input$seurat_cluster_set.dd,input$gene_of_interest.dd)) %>%
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

    if(FALSE)
      input <- list(features_of_interest.tb={Matrix::rowSums(seurat) %>% (function(x) x[x>100]) %>% names()}, seurat_cluster_set.dd='integrated_snn_res.1.2')
    
    progress$inc(detail='Getting list of feature IDs')
    feature_names <- input$features_of_interest.tb
    if(feature_names[1]=='')
      Matrix::rowSums(seurat) %>%
        sort() %>%
        head(2000) %>%
        sample(size=200) %>%
        names() -> feature_names
    feature_names %<>% str_split(pattern='\\s|;|,') %>% unlist() %>% (function(x) x[x %in% rownames(seurat)])

    progress$inc(detail='Collecting feature values matrix')
    FetchData(object=seurat, vars=feature_names) %>%
      rownames_to_column('cell_id') -> feature_values_matrix
    FetchData(object=seurat, vars=input$features_heatmap.seurat_cluster_set.dd) %>%
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
                           width=unit(1, 'npc'),
                           height=unit(1, 'npc')))}) %>% debounce(2500)
  
  # ###############################################################################################
  # return plots to shiny -------------------------------------------------------------------------
  progress$inc(detail='Returning plots to Shiny')

  ## highlighting genes tab
  renderPlot(genes_highlighting.clustered_map.plot()) -> output$`genes_highlighting-clustered_map`
  renderPlot(genes_highlighting.gene_expression_map.plot()) -> output$`genes_highlighting-gene_expression_map`
  renderPlot(genes_highlighting.expression_per_cluster.plot()) -> output$`genes_highlighting-expression_per_cluster`

  output$genes_highlighting.n_cells_box <- renderValueBox(expr={
    valueBox(value=scales::comma(ncol(seurat)),
             subtitle='Cells in map',
             icon=icon('galactic-republic'),
             color='purple')})
  output$genes_highlighting.n_clusters_box <- renderValueBox(expr={
    valueBox(value=scales::comma(seurat_cluster_per_set[input$seurat_cluster_set.dd]),
             subtitle='Cell clusters',
             icon=icon('first-order'),
             color='purple')})
  output$genes_highlighting.selected_gene_box <- renderValueBox(expr={
    valueBox(value=input$gene_of_interest.dd,
             subtitle={gene_names_to_description[input$gene_of_interest.dd] %>% str_trunc(width=100) %>% (function(x) if_else(is.na(x), 'Selected gene', x))},
             icon=icon('jedi-order'),
             color='purple')})
  output$genes_highlighting.n_genes_box <- renderValueBox(expr={
    valueBox(value=scales::comma(sapply(seurat@assays, function(x) nrow(x@data)) %>% max()),
             subtitle='Unique genes detected (max)',
             icon=icon('galactic-senate'),
             color='purple')})
  output$genes_highlighting.n_reads_box <- renderValueBox(expr={
    valueBox(value=scales::comma(sum(seurat$nCount_RNA)),
             subtitle='Total reads in cells',
             icon=icon('old-republic'),
             color='purple')})
  output$genes_highlighting.n_reads_per_cell_box <- renderValueBox(expr={
    valueBox(value=scales::comma(round(median(seurat$nCount_RNA), digits=1)),
             subtitle='Median reads per cell',
             icon=icon('frog'),
             color='purple')})
  output$genes_highlighting.n_genes_per_cell_box <- renderValueBox(expr={
    valueBox(value=scales::comma(round(median(seurat$nFeature_RNA), digits=1)),
             subtitle='Median genes per cell',
             icon=icon('crow'),
             color='purple')})

  ## cell filtering tab
  renderPlot(cell_filtering.total_expression_knee.plot()) -> output$`cell_filtering-total_expression_knee`
  renderPlot(cell_filtering.unique_genes_knee.plot()) -> output$`cell_filtering-unique_genes_knee`
  renderPlot(cell_filtering.percent_mitochondria_knee.plot()) -> output$`cell_filtering-percent_mitochondria_knee`

  renderPlot(cell_filtering.total_expression_density.plot()) -> output$`cell_filtering-total_expression_density`
  renderPlot(cell_filtering.unique_features_density.plot()) -> output$`cell_filtering-unique_features_density`
  renderPlot(cell_filtering.percent_mitochondria_density.plot()) -> output$`cell_filtering-percent_mitochondria_density`

  output$cell_filtering.n_reads_box <- renderValueBox({
    react_to_cell_filtering()
    valueBox(value=scales::comma(cell_filtering_data.reactions$total_reads),
             subtitle=sprintf(fmt='Total reads remaining (%.1f%%)', cell_filtering_data.reactions$total_reads/cell_filtering_data.reference$total_reads*100),
             icon=icon('old-republic'),
             color='purple')})
  output$cell_filtering.n_cells_box <- renderValueBox({
    react_to_cell_filtering()
    valueBox(value=scales::comma(cell_filtering_data.reactions$n_cells),
             subtitle=sprintf(fmt='Cells remaining (%.1f%%)', cell_filtering_data.reactions$n_cells/cell_filtering_data.reference$n_cells*100),
             icon=icon('galactic-republic'),
             color='purple')})
  output$cell_filtering.n_reads_per_cell_box <- renderValueBox({
    react_to_cell_filtering()
    valueBox(value=scales::comma(cell_filtering_data.reactions$median_reads_per_cell),
             subtitle=sprintf(fmt='Median reads per cell (%+d)', cell_filtering_data.reactions$median_reads_per_cell-cell_filtering_data.reference$median_reads_per_cell),
             icon=icon('frog'),
             color='purple')})
  output$cell_filtering.n_genes_per_cell_box <- renderValueBox({
    react_to_cell_filtering()
    valueBox(value=scales::comma(cell_filtering_data.reactions$median_genes_per_cell),
             subtitle=sprintf(fmt='Median genes per cell (%+d)', cell_filtering_data.reactions$median_genes_per_cell-cell_filtering_data.reference$median_genes_per_cell),
             icon=icon('crow'),
             color='purple')})

  # features heatmap tab
  renderPlot(features_heatmap.heatmap.plot()) -> output$`features_heatmap-heatmap`
  
  # sidebar
  ## dynamic sidebar outputs can be listed here
  
  # shared text boxes
  
  ## project name
  project_name_box_opts <- list(value={seurat@project.name %>% str_replace_all('_', ' ') %>% str_to_upper()}, subtitle='Loaded Seurat object', icon=icon('certificate'), color='purple')
  output$cell_filtering.project_name_box <- renderValueBox(expr={do.call(what=valueBox, args=project_name_box_opts)})
  output$genes_highlighting.project_name_box <- renderValueBox(expr={do.call(what=valueBox, args=project_name_box_opts)})
  output$features_heatmap.project_name_box <- renderValueBox(expr={do.call(what=valueBox, args=project_name_box_opts)})
}
