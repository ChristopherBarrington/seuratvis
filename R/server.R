
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

  ## call reduced dimensionality plot modules
  for(id in module_environments$reduced_dimension_plots$id)
    callModule(module=reduced_dimension_plot.server, id=id)

  ## call feature values per cluster plot modules
  for(id in module_environments$feature_values_per_cluster_plots$id)
    callModule(module=feature_values_per_cluster_plot.server, id=id)

  # ###############################################################################################
  # return plots to shiny -------------------------------------------------------------------------
  progress$inc(detail='Returning plots to Shiny')

  ## highlighting genes tab

  callModule(module=project_name_text_box.server, id='gene_highlighting')
  callModule(module=number_of_cells_text_box.server, id='gene_highlighting')
  callModule(module=number_of_clusters_text_box.server, id='gene_highlighting')
  callModule(module=gene_name_and_description_text_box.server, id='gene_highlighting')
  callModule(module=number_of_reads_text_box.server, id='gene_highlighting')
  callModule(module=number_of_genes_in_assay_text_box.server, id='gene_highlighting')
  callModule(module=number_of_reads_per_cell_text_box.server, id='gene_highlighting')
  callModule(module=number_of_genes_per_cell_text_box.server, id='gene_highlighting')

  ## cell filtering tab

  callModule(module=project_name_text_box.server, id='cell_filtering')
  callModule(module=number_of_reads_text_box.server, id='cell_filtering')
  callModule(module=number_of_cells_text_box.server, id='cell_filtering')
  callModule(module=number_of_reads_per_cell_text_box.server, id='cell_filtering')
  callModule(module=number_of_genes_per_cell_text_box.server, id='cell_filtering')

  # sidebar
  ## dynamic sidebar outputs can be listed here

  # any code to exectue when the session ends
  session$onSessionEnded(function() {message('### session ended')})
}
