#'
#'
dataset_filtering.server <- function(input, output, session, seurat, filters=list()) {
  reactiveValues() -> filtering_parameters

  # update the reactive when the filters are changed
  observe(label='dataset_filtering/observe', priority=-999, x={
    req(seurat$metadata)

    # if there are no filters provided
    if(length(filters)==0)
      return(NULL)

    # if any of the filters are not uet initialised
    if(lapply(filters, reactiveValuesToList) %>% sapply(length) %>% equals(0) %>% any())
      return(NULL)

    # if any of the variables are null
    if(lapply(filters, reactiveValuesToList) %>% sapply(function(x) is.null(pluck(x, 'variable'))) %>% any())
      return(NULL)

    # if any of the variables are missing from the metadata
    if(!lapply(filters, reactiveValuesToList) %>% sapply(function(x) pluck(x, 'variable')) %>% is.element(set=names(seurat$metadata)) %>% all())
      return(NULL)

    # get the values in the list of reactives

    lapply(filters, reactiveValuesToList)  %>%
      lapply(function(x) x %>% extract(str_detect(string=names(x), pattern='variable|min|max|in_set'))) %>% # pick out the elements we can use
      plyr::ldply(as.data.frame) %>%
      gather(key=logic, value=value, -variable) %>%
      drop_na() -> filters_df

    # prepare a condition to filter the metadata
    ## by default, select all cells (no filtering)
    filter_condition <- 'TRUE'

    ## if there are some filters, prepare a condition to filter the cells
    if(nrow(filters_df)>0)
      filters_df %>%
        # gather(key=logic, value=value, -variable) %>%
        # drop_na() %>%
        mutate(logic=factor(logic, levels=c('min', 'max', 'in_set')),
               logic=fct_recode(logic, `>=`='min', `<=`='max', ` %in% `='in_set'),
               value=str_trim((value))) %>%
        apply(1, str_c, collapse='') %>%
        str_c(collapse=' & ') -> filter_condition

    # filter the cells
    filtered_cell_metadata <- filter(.data=seurat$metadata, eval(parse(text=filter_condition)))

    # pull out the filtered variables
    filtered_n_umi_values <- pluck(filtered_cell_metadata, seurat$n_umi_variable)
    filtered_n_features_values <- pluck(filtered_cell_metadata, seurat$n_features_variable)
    filtered_proportion_mt_values <- pluck(filtered_cell_metadata, seurat$proportion_mt_variable)

    # save the filtered data to the reactive
    filtering_parameters$n_cells <- nrow(filtered_cell_metadata)
    filtering_parameters$n_umi_values <- filtered_n_umi_values
    filtering_parameters$n_umi <- sum(filtered_n_umi_values)
    filtering_parameters$n_features_values <- filtered_n_features_values
    filtering_parameters$done <- rnorm(1)})

  # return the reactive values list
  return(filtering_parameters)
}
