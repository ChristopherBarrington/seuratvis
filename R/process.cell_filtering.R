#' React to filtering thresholds and filter Seurat
#' 
#' @details
#' Updates the Seurat object and thresholds reactive values.
#' 
#' @rdname cell_filtering
#'
cell_filtering.server <- function(input, output, session, seurat, ...) {
  message('### cell_filtering.server')

  # get environments containing variables to run/configure this object
  collect_environments(id=parent.frame()$id, module='total_umi_per_cell_filter') # provides `seuratvis_env`, `server_env` and `module_env`
  session_server <- get(x='session', env=server_env)

  filtering.rv <- reactiveValues()

  # react to the filtering parameters being changed
  observeEvent(eventExpr=filtering.rv$updated_parameter, handlerExpr={
    # make sure these elements are defined
    req(seurat_object.reactions$seurat)
    req(seurat$object)

    req(seurat$n_features_values)
    req(seurat$n_umi_values)
    req(seurat$proportion_mt_values)

    req(seurat$n_features_variable)
    req(seurat$n_umi_variable)

    # create variables for shorthand
    cell_metadata <- seurat$metadata
    min_umi_per_cell <- filtering.rv$total_umi_per_cell_min
    max_umi_per_cell <- filtering.rv$total_umi_per_cell_max
    min_features_per_cell <- filtering.rv$features_per_cell_min
    max_features_per_cell <- filtering.rv$features_per_cell_max
    max_percent_mitochondria <- filtering.rv$max_percent_mitochondria

    if(!all(min_umi_per_cell>0 & max_umi_per_cell>0 &
            min_features_per_cell>0 & max_features_per_cell>0 &
            max_percent_mitochondria>0))
      return(NULL)

    # send a message
    message('### cell_filtering.server-observeEvent-filtering.rv$updated_parameter')

    n_features_values <- seurat$n_features_values %>% unlist(use.names=FALSE)
    proportion_mt_values <- seurat$proportion_mt_values %>% unlist(use.names=FALSE)
    n_umi_values <- seurat$n_umi_values %>% unlist(use.names=FALSE)

    # filter Seurat object
    ({n_features_values %>% between(left=min_features_per_cell, right=max_features_per_cell)} &
     {n_umi_values %>% between(left=min_umi_per_cell, right=max_umi_per_cell)} &
     {proportion_mt_values<=max_percent_mitochondria}) %>%
      filter(.data=cell_metadata) -> filtered_cell_metadata

    filtered_n_features_values <- pluck(filtered_cell_metadata, seurat$n_features_variable)
    filtered_n_umi_values <- pluck(filtered_cell_metadata, seurat$n_umi_variable)

    # save values to filtering reactive
    filtering.rv$n_cells <- nrow(filtered_cell_metadata)
    filtering.rv$n_umi <- sum(filtered_n_umi_values)
    filtering.rv$n_features_values <- filtered_n_features_values
    filtering.rv$n_umi_values <- filtered_n_umi_values
    filtering.rv$done <- rnorm(1)})

  # react to a new seurat being loaded or variable names being changed
  ## proportion mitochondria
  # observeEvent(eventExpr=seurat$proportion_mt_values_max, handlerExpr={filtering.rv$max_percent_mitochondria <- seurat$proportion_mt_values_max})

  ## number of features
  # observeEvent(eventExpr=seurat$n_features_values_min, handlerExpr={filtering.rv$features_per_cell_min <- seurat$n_features_values_min})
  # observeEvent(eventExpr=seurat$n_features_values_max, handlerExpr={filtering.rv$features_per_cell_max <- seurat$n_features_values_max})
  observeEvent(eventExpr=seurat$n_features_values, handlerExpr={filtering.rv$n_features_values <- seurat$n_features_values})

  ## number of cells
  observeEvent(eventExpr=seurat$n_cells, handlerExpr={filtering.rv$n_cells <- seurat$n_cells})

  ## number of umi
  # observeEvent(eventExpr=seurat$n_umi_values_min, handlerExpr={filtering.rv$total_umi_per_cell_min <- seurat$n_umi_values_min})
  # observeEvent(eventExpr=seurat$n_umi_values_max, handlerExpr={filtering.rv$total_umi_per_cell_max <- seurat$n_umi_values_max})
  observeEvent(eventExpr=seurat$n_umi_values, handlerExpr={filtering.rv$n_umi_values <- seurat$n_umi_values})
  observeEvent(eventExpr=seurat$n_umi, handlerExpr={filtering.rv$n_umi <- seurat$n_umi})

  # return the reactive values list
  return(filtering.rv)
}
