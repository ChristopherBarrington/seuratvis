#' 
#' 
dataset_info_text_box.ui <- function(id, width=12)
  valueBoxOutput(outputId=NS(id, 'box'), width=width)

#'
#' @details Requires font-awesome version 5+ to be installed in shiny/www/shared/fontawesome
#' 
dataset_info_text_box.defaults <- function(...)
  list(color='purple',
       icon=sample(x=c('paw','otter','hippo','dog','spider','kiwi-bird','horse-head','horse','frog','fish','dragon','dove','crow','cat'), size=1) %>% icon())

#'
#'  
dataset_info_text_box.project_name <- function(input, output, session, seurat) {
  renderValueBox(expr={
    req(seurat$formatted_project_name)

    list(value=seurat$formatted_project_name,
         subtitle='Loaded Seurat object') %>%
      modifyList(x=dataset_info_text_box.defaults()) %>%
      do.call(what=valueBox)}) -> output$box
}

#'
#' @import scales
#'  
dataset_info_text_box.n_umi <- function(input, output, session, seurat) {
  renderValueBox(expr={
    req(seurat$n_umi_sum)

    list(value=comma(seurat$n_umi_sum),
         subtitle='Total UMI in cells') %>%
      modifyList(x=dataset_info_text_box.defaults()) %>%
      do.call(what=valueBox)}) -> output$box
}

#'
#' @import scales
#'  
dataset_info_text_box.n_filtered_umi <- function(input, output, session, seurat, cell_filtering) {
  renderValueBox(expr={
    req(seurat$n_umi_sum)
    req(cell_filtering$n_umi)

    n_reference <- seurat$n_umi_sum
    n_filtered <- cell_filtering$n_umi
    n_removed <- n_reference-n_filtered
    subtitle <- sprintf(fmt='%s UMI removed (%.1f%% remain)', comma(n_removed), n_filtered/n_reference*100)

    list(value=comma(n_filtered),
         subtitle=subtitle) %>%
      modifyList(x=dataset_info_text_box.defaults()) %>%
      do.call(what=valueBox)}) -> output$box
}

#'
#' @import scales
#'  
dataset_info_text_box.n_cells <- function(input, output, session, seurat) {
  renderValueBox(expr={
    req(seurat$n_cells)

    list(value=comma(seurat$n_cells),
         subtitle='Cells in map') %>%
      modifyList(x=dataset_info_text_box.defaults()) %>%
      do.call(what=valueBox)}) -> output$box
}

#'
#' @import scales
#'  
dataset_info_text_box.n_filtered_cells <- function(input, output, session, seurat, cell_filtering) {
  renderValueBox(expr={
    req(seurat$n_cells)
    req(cell_filtering$n_cells)

    n_reference <- seurat$n_cells
    n_filtered <- cell_filtering$n_cells
    n_removed <- n_reference-n_filtered
    subtitle <- sprintf(fmt='%s cells removed (%.1f%% remain)', comma(n_removed), n_filtered/n_reference*100)

    list(value=n_filtered %>% scales::comma(),
         subtitle=subtitle) %>%
      modifyList(x=dataset_info_text_box.defaults()) %>%
      do.call(what=valueBox)}) -> output$box
}

#'
#' @import scales
#'  
dataset_info_text_box.n_umi_per_cell <- function(input, output, session, seurat, sf=median, sf_name='Median') {
  renderValueBox(expr={
    req(seurat$n_umi_values)

    n_reference <- seurat$n_umi_values %>% unlist() %>% sf()
    subtitle <- sprintf('%s reads per cell', str_to_title(sf_name))

    list(value=n_reference %>% scales::comma(), subtitle=subtitle) %>%
      modifyList(x=dataset_info_text_box.defaults()) %>%
      do.call(what=valueBox)}) -> output$box
}

#'
#' @import scales
#'  
dataset_info_text_box.n_umi_per_filtered_cell <- function(input, output, session, seurat, cell_filtering, sf=median, sf_name='Median') {
  renderValueBox(expr={
    req(seurat$n_umi_values)
    req(cell_filtering$n_umi_values)

    n_reference <- seurat$n_umi_values %>% unlist() %>% sf()
    n_filtered <- cell_filtering$n_umi_values %>% unlist() %>%  (function(x) ifelse(is.null(x), 0, sf(x)))
    subtitle <- sprintf(fmt='%s reads per cell (%s)', str_to_title(sf_name), comma(n_filtered-n_reference) %>% ifelse(str_detect(., '^-'), ., str_c('+', .)))

    list(value=n_filtered %>% scales::comma(), subtitle=subtitle) %>%
      modifyList(x=dataset_info_text_box.defaults()) %>%
      do.call(what=valueBox)}) -> output$box
}

#'
#' @import scales
#'  
dataset_info_text_box.n_features_per_cell <- function(input, output, session, seurat, sf=median, sf_name='Median') {
  renderValueBox(expr={
    req(seurat$n_features_values)

    n_reference <- seurat$n_features_values %>% unlist() %>% sf()
    subtitle <- sprintf('%s features per cell', str_to_title(sf_name))

    list(value=n_reference %>% scales::comma(), subtitle=subtitle) %>%
      modifyList(x=dataset_info_text_box.defaults()) %>%
      do.call(what=valueBox)}) -> output$box
}

#'
#' @import scales
#'  
dataset_info_text_box.n_features_per_filtered_cell <- function(input, output, session, seurat, cell_filtering, sf=median, sf_name='Median') {
  renderValueBox(expr={
    req(seurat$n_features_values)
    req(cell_filtering$n_features_values)

    n_reference <- seurat$n_features_values %>% unlist() %>% sf()
    n_filtered <- cell_filtering$n_features_values %>% unlist() %>%  (function(x) ifelse(is.null(x), 0, sf(x)))
    subtitle <- sprintf(fmt='%s features per cell (%s)', str_to_title(sf_name), comma(n_filtered-n_reference) %>% ifelse(str_detect(., '^-'), ., str_c('+', .)))

    list(value=n_filtered %>% scales::comma(), subtitle=subtitle) %>%
      modifyList(x=dataset_info_text_box.defaults()) %>%
      do.call(what=valueBox)}) -> output$box
}

