#' 
#' 
filter_proportion_mt.ui <- function(id, label='Proportion Mt', seurat) {
  value <- seurat$proportion_mt_values_max %>% add(0.05) %>% round(digits=1)

  # make ui elements
  div(tags$h6('Maximum', style='display: inline;'),
      numericInput(inputId=NS(id, 'max_percent_mt'), label=NULL, value=value, min=0, max=value, step=0.1, width='100%')) -> high_ui

  # return ui element(s)
  tagList(tags$label(label), br(), high_ui)
}

#'
#'
filter_proportion_mt.server <- function(input, output, session, cell_filtering) {
  # react to the input elements
  ## max proportion mitochondrial umi
  observeEvent(eventExpr=input$max_percent_mt, ignoreInit=TRUE, handlerExpr={
    # update the reactive
    value <- input$max_percent_mt
    if(value != sprintf(fmt='%.1f', value)) # if the value is not already 1dp formatted, reformat it
      value %<>% add(0.05) %>% round(digits=1)

    cell_filtering$proportion_mt_max <- value
    cell_filtering$updated_parameter <- rnorm(1)})

  # react to the cell filtering parameter being changed eg by the brush
  ## max proportion mitochondrial umi
  observeEvent(eventExpr=cell_filtering$proportion_mt_max, handlerExpr={
    req(input$max_percent_mt)
    high <- cell_filtering$proportion_mt_max
    if(high != sprintf(fmt='%.1f', high)) # if the value is not already 1dp formatted, reformat it
      high <- cell_filtering$proportion_mt_max %>% add(0.05) %>% round(digits=1)
    if(high != input$max_percent_mt)
      updateNumericInput(session=session, inputId='max_percent_mt', value=high)})
}
