
10^(0:9) -> major_breaks_log10
(2:9) * rep(major_breaks_log10, each=8) -> minor_breaks_log10

module_environments <- new.env()
seurat_object.reactions <- reactiveValues()
filtering_parameters.reactions <- reactiveValues()
filtered_cells.reactions <- reactiveValues()

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

  ## react to opening tab with a filtered object loaded
  observeEvent(input$sidebarmenu, {
    if(!is.null(seurat_object.reactions$seurat) & input$sidebarmenu=='cell_filtering-tab' && (!is.null(seurat_object.reactions$seurat@misc$cells_filtered) && seurat_object.reactions$seurat@misc$cells_filtered))
      sendSweetAlert(session=session, type='success', html=TRUE,
                     title='Notice', btn_labels='Great!',
                     text=tags$span('It looks like low-quality cells have already been removed from this Seurat object:', tags$h5(tags$code('@misc$cells_filtered == TRUE'))),
                     closeOnClickOutside=TRUE, showCloseButton=FALSE)})

  ## load the filter_seurat module
  callModule(module=cell_filtering.server, id='seuratvis')

  ## call knee plot modules
  for(id in module_environments$knee_plots$id)
    callModule(module=knee_plot.server, id=id)

  ## call boxplot modules
  for(id in module_environments$boxplot_plots$id)
    callModule(module=boxplot_plot.server, id=id)

  ## call density plot modules
  for(id in module_environments$density_plots$id)
    callModule(module=density_plot.server, id=id)

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
      theme(legend.position='none') + labs(title=seurat_object.reactions$picked_feature)
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
      theme(legend.position='bottom', legend.title=element_blank()) + labs(title=seurat_object.reactions$picked_feature)
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

  # features heatmap tab
  callModule(project_name_text_box.server, id='features_heatmap')
  renderPlot(features_heatmap.heatmap.plot()) -> output$`features_heatmap-heatmap`

  # sidebar
  ## dynamic sidebar outputs can be listed here

  # any code to exectue when the session ends
  session$onSessionEnded(function() {message('### session ended')})
}
