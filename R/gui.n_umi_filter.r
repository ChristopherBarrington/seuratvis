#' 
#' 
filter_n_umi.ui <- function(id, label='UMI per cell', seurat, low=TRUE, high=TRUE) {
  # make ui elements
  div(tags$h6('Minimum', style='display: inline;'),
      numericInput(inputId=NS(id, 'min_umi'), label=NULL, value=seurat$n_umi_values_min, min=seurat$n_umi_values_min, step=50, width='100%')) -> low_ui

  div(tags$h6('Maximum', style='display: inline;'),
      numericInput(inputId=NS(id, 'max_umi'), label=NULL, value=seurat$n_umi_values_max, max=seurat$n_umi_values_max, step=50, width='100%')) -> high_ui

  # return ui element(s)
  tagList(tags$label(label), br(), if(low) low_ui, if(high) high_ui)
}

#'
#'
filter_n_umi.server <- function(input, output, session, seurat, cell_filtering, linked_density_plot) {
  filtering <- reactiveValues()

  # react to the input elements
  ## min_umi
  observeEvent(eventExpr=input$min_umi, handlerExpr={
    # check that max>min and update the reactive
    if(input$min_umi<=input$max_umi) {
      cell_filtering$n_umi_min <- round(input$min_umi, digits=0)
      filtering$min <- round(input$min_umi, digits=0)
    } else {
      updateNumericInput(session=session, inputId='min_umis', value=input$max_umis)
    }})

  ## max_umi
  observeEvent(eventExpr=input$max_umi, handlerExpr={
    # check that max>min and update the reactive
    if(input$max_umi>=input$min_umi) {
      cell_filtering$n_umi_max <- round(input$max_umi, digits=0)
      filtering$max <- round(input$max_umi, digits=0)
    } else {
      updateNumericInput(session=session, inputId='max_umis', value=input$min_umis)
    }})

  # react to the cell filtering parameter being changed eg by the brush
  ## min umi
  observeEvent(eventExpr=linked_density_plot$min, ignoreInit=TRUE, handlerExpr={
    updateNumericInput(session=session, inputId='min_umi', value=linked_density_plot$min)})

  ## max umi
  observeEvent(eventExpr=linked_density_plot$max, ignoreInit=TRUE, handlerExpr={
    updateNumericInput(session=session, inputId='max_umi', value=linked_density_plot$max)})

  # react to a new seurat being loaded
  observeEvent(eventExpr=seurat$n_umi_updated, label='filter_n_umi/object', handlerExpr={
    filtering$variable <- seurat$n_umi_variable
    filtering$min <- seurat$n_umi_values_min
    filtering$max <- seurat$n_umi_values_max})

  # return the reactiveValues
  return(filtering)
}
