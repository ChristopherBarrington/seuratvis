#' React to filtering thresholds and filter Seurat
#' 
#' @details
#' Updates the Seurat object and thresholds reactive values.
#' 
#' @rdname cell_filtering
#'
cell_filtering.server <- function(input, output, session) {
  message('### cell_filtering.server')

  # get environments containing variables to run/configure this object
  collect_environments(id=parent.frame()$id, module='total_umi_per_cell_filter') # provides `seuratvis_env`, `server_env` and `module_env`
  session_server <- get(x='session', env=server_env)

  # react to the filtering parameters being changed
  observeEvent(eventExpr=reactiveValuesToList(filtering_parameters.reactions), handlerExpr={
    message('### cell_filtering.server-observeEvent-reactiveValuesToList(filtered_cells.reactions)')

    # make sure seurat object is loaded
    req(seurat_object.reactions$seurat)

    # create variables for shorthand
    cell_metadata <- seurat_object.reactions$cell_metadata
    min_umi_per_cell <- filtering_parameters.reactions$total_umi_per_cell_min
    max_umi_per_cell <- filtering_parameters.reactions$total_umi_per_cell_max
    min_features_per_cell <- filtering_parameters.reactions$features_per_cell_min
    max_features_per_cell <- filtering_parameters.reactions$features_per_cell_max
    max_percent_mitochondria <- filtering_parameters.reactions$max_percent_mitochondria

    # filter Seurat object
    cell_metadata %>%
      filter(percent_mt <= max_percent_mitochondria &
             between(x=nFeature_RNA, left=min_features_per_cell, right=max_features_per_cell) &
             between(x=nCount_RNA, left=min_umi_per_cell, right=max_umi_per_cell)) -> filtered_cell_metadata

    # save values to filtering reactive
    filtered_cells.reactions$n_cells <- nrow(filtered_cell_metadata)
    filtered_cells.reactions$n_umi <- sum(filtered_cell_metadata$nCount_RNA)
    filtered_cells.reactions$median_umi_per_cell <- round(x=median(filtered_cell_metadata$nCount_RNA), digits=0)
    filtered_cells.reactions$median_features_per_cell <- round(x=median(filtered_cell_metadata$nFeature_RNA), digits=0)
    filtered_cells.reactions$cell_metadata <- filtered_cell_metadata})

  # initialise filtering thresholds and filtered object reactives when Seurat object is loaded
  observeEvent(eventExpr=seurat_object.reactions$seurat, handlerExpr={
    sprintf(fmt='### cell_filtering.server-observeEvent-seurat_object.reactions$seurat [%s]', seurat_object.reactions$formatted.project.name) %>% message()

    # create variables for shorthand
    seurat <- seurat_object.reactions$seurat
    cell_metadata <- seurat_object.reactions$cell_metadata

    # get the initialisation values before filtering
    list(project=Project(seurat),
         n_cells=nrow(cell_metadata),
         n_umi=sum(cell_metadata$nCount_RNA),
         
         total_umi_per_cell_min=min(cell_metadata$nCount_RNA),
         total_umi_per_cell_max=max(cell_metadata$nCount_RNA),
         
         features_per_cell_min=min(cell_metadata$nFeature_RNA),
         features_per_cell_max=max(cell_metadata$nFeature_RNA),
         
         max_percent_mitochondria=round(max(cell_metadata$percent_mt)+0.05, digits=1)) -> initial_values

    # save the values into the reactives
    for(i in names(initial_values))
      filtering_parameters.reactions[[i]] <- initial_values[[i]]
    filtered_cells.reactions$cell_metadata <- cell_metadata})
}
