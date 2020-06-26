#'
#'
dataset_filtering.server <- function(input, output, session, seurat) {
  reactiveValues() -> filtering_parameters

  # update the reactive when the filters are changed
  observe(label='dataset_filtering/reaction', x={
    req(filtering_parameters$metadata)
    req(filtering_parameters$n_umi_min)
    req(filtering_parameters$n_umi_max)
    req(filtering_parameters$n_features_min)
    req(filtering_parameters$n_features_max)
    req(filtering_parameters$proportion_mt_max)

    req(seurat$n_umi_values)
    req(seurat$n_features_values)
    req(seurat$proportion_mt_values)
    req(seurat$n_umi_variable)
    req(seurat$n_features_variable)
    req(seurat$proportion_mt_variable)

    # copy the variables
    cell_metadata <- filtering_parameters$metadata
    min_umi_per_cell <- filtering_parameters$n_umi_min
    max_umi_per_cell <- filtering_parameters$n_umi_max
    min_features_per_cell <- filtering_parameters$n_features_min
    max_features_per_cell <- filtering_parameters$n_features_max
    max_percent_mitochondria <- filtering_parameters$proportion_mt_max

    # check that all thresholds are non-zero
    if(!all(min_umi_per_cell>0 & max_umi_per_cell>0 &
            min_features_per_cell>0 & max_features_per_cell>0 &
            max_percent_mitochondria>0))
      return(NULL)

    # filter the value vectors
    ({seurat$n_umi_values %>% between(left=min_umi_per_cell, right=max_umi_per_cell)} &
     {seurat$n_features_values %>% between(left=min_features_per_cell, right=max_features_per_cell)} &
     {seurat$proportion_mt_values<=max_percent_mitochondria}) %>%
      filter(.data=cell_metadata) -> filtered_cell_metadata

    # pull out the filtered variables
    filtered_n_umi_values <- pluck(filtered_cell_metadata, seurat$n_umi_variable)
    filtered_n_features_values <- pluck(filtered_cell_metadata, seurat$n_features_variable)
    filtered_proportion_mt_values <- pluck(filtered_cell_metadata, seurat$proportion_mt_variable)

    # save the filtered data to the reactive
    filtering_parameters$n_cells <- nrow(filtered_cell_metadata)
    filtering_parameters$n_umi <- sum(filtered_n_umi_values)
    filtering_parameters$n_features_values <- filtered_n_features_values
    filtering_parameters$n_umi_values <- filtered_n_umi_values
    filtering_parameters$done <- rnorm(1)})

  # react to a new seurat being loaded
  observeEvent(eventExpr=c(seurat$n_features_updated, seurat$n_umi_updated, seurat$proportion_mt_updated), label='dataset_filtering/object', handlerExpr={
    req(seurat$n_features_updated)
    req(seurat$n_umi_updated)
    req(seurat$proportion_mt_updated)

    # wipe the reactive values
    for(i in names(filtering_parameters)) filtering_parameters[[i]] <- NULL

    # set the starting thresholds
    filtering_parameters$n_umi_min <- seurat$n_umi_values_min
    filtering_parameters$n_umi_max <- seurat$n_umi_values_max
    filtering_parameters$n_features_min <- seurat$n_features_values_min
    filtering_parameters$n_features_max <- seurat$n_features_values_max
    filtering_parameters$proportion_mt_max <- seurat$proportion_mt_values_max})

  # return the reactive values list
  return(filtering_parameters)
}
