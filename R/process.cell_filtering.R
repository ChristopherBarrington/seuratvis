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
    # make sure these elements are defined
    req(seurat_object.reactions$seurat)

    req(seurat_object.reactions$n_features_values)
    req(seurat_object.reactions$n_umi_values)
    req(seurat_object.reactions$proportion_mt_values)

    req(seurat_configuration.reactions$n_features_variable)
    req(seurat_configuration.reactions$n_umi_variable)

    # send a message
    message('### cell_filtering.server-observeEvent-reactiveValuesToList(filtered_cells.reactions)')

    # create variables for shorthand
    cell_metadata <- seurat_object.reactions$cell_metadata
    min_umi_per_cell <- filtering_parameters.reactions$total_umi_per_cell_min
    max_umi_per_cell <- filtering_parameters.reactions$total_umi_per_cell_max
    min_features_per_cell <- filtering_parameters.reactions$features_per_cell_min
    max_features_per_cell <- filtering_parameters.reactions$features_per_cell_max
    max_percent_mitochondria <- filtering_parameters.reactions$max_percent_mitochondria

    if(!all(min_umi_per_cell>0 & max_umi_per_cell>0 &
            min_features_per_cell>0 & max_features_per_cell>0 &
            max_percent_mitochondria>0))
      return(NULL)

    n_features_values <- seurat_object.reactions$n_features_values %>% unlist(use.names=FALSE)
    proportion_mt_values <- seurat_object.reactions$proportion_mt_values %>% unlist(use.names=FALSE)
    n_umi_values <- seurat_object.reactions$n_umi_values %>% unlist(use.names=FALSE)

    # filter Seurat object
    ({n_features_values %>% between(left=min_features_per_cell, right=max_features_per_cell)} &
     {n_umi_values %>% between(left=min_umi_per_cell, right=max_umi_per_cell)} &
     {proportion_mt_values<=max_percent_mitochondria}) %>%
      filter(.data=cell_metadata) -> filtered_cell_metadata

    filtered_n_features_values <- pluck(filtered_cell_metadata, seurat_configuration.reactions$n_features_variable)
    filtered_n_umi_values <- pluck(filtered_cell_metadata, seurat_configuration.reactions$n_umi_variable)

    # save values to filtering reactive
    filtered_cells.reactions$n_cells <- nrow(filtered_cell_metadata)
    filtered_cells.reactions$n_umi <- sum(filtered_n_umi_values)
    filtered_cells.reactions$n_features_values <- filtered_n_features_values
    filtered_cells.reactions$n_umi_values <- filtered_n_umi_values
    filtered_cells.reactions$median_features_per_cell <- round(x=median(filtered_n_features_values), digits=0)
    filtered_cells.reactions$median_umi_per_cell <- round(x=median(filtered_n_umi_values), digits=0)
    filtered_cells.reactions$cell_metadata <- filtered_cell_metadata})
}
