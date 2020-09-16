#' 
#' 
filter_n_features.ui <- function(id, label='Features per cell', seurat, low=TRUE, high=TRUE) {
  # make ui elements
  div(tags$h6('Minimum', style='display: inline;'),
      numericInput(inputId=NS(id, 'min_features'), label=NULL, value=seurat$n_features_values_min, min=seurat$n_features_values_min, step=50, width='100%')) -> low_ui

  div(tags$h6('Maximum', style='display: inline;'),
      numericInput(inputId=NS(id, 'max_features'), label=NULL, value=seurat$n_features_values_max, max=seurat$n_features_values_max, step=50, width='100%')) -> high_ui

  # return ui element(s)
  tagList(tags$label(label), br(), if(low) low_ui, if(high) high_ui)
}

#'
#'
filter_n_features.server <- function(input, output, session, seurat, cell_filtering, linked_density_plot) {
  filtering <- reactiveValues()

  # react to the input elements
  ## min features
  observeEvent(eventExpr=input$min_features, handlerExpr={
    # check that max>min and update the reactive
    if(input$min_features<=input$max_features) {
      cell_filtering$n_features_min <- round(input$min_features, digits=0)
      filtering$min <- round(input$min_features, digits=0)
    } else {
      updateNumericInput(session=session, inputId='min_features', value=input$max_features)
    }})

  ## max features
  observeEvent(eventExpr=input$max_features, handlerExpr={
    # check that max>min and update the reactive
    if(input$max_features>=input$min_features) {
      cell_filtering$n_features_max <- round(input$max_features, digits=0)
      filtering$max <- round(input$max_features, digits=0)
    } else {
      updateNumericInput(session=session, inputId='max_features', value=input$min_features)
    }})

  # react to the cell filtering parameter being changed eg by the brush
  ## min features
  observeEvent(eventExpr=linked_density_plot$min, ignoreInit=TRUE, handlerExpr={
    updateNumericInput(session=session, inputId='min_features', value=linked_density_plot$min)})

  ## max features
  observeEvent(eventExpr=linked_density_plot$max, ignoreInit=TRUE, handlerExpr={
    updateNumericInput(session=session, inputId='max_features', value=linked_density_plot$max)})

  # react to a new seurat being loaded
  observeEvent(eventExpr=seurat$n_features_updated, label='filter_n_features/object', handlerExpr={
    filtering$variable <- seurat$n_features_variable
    filtering$min <- seurat$n_features_values_min
    filtering$max <- seurat$n_features_values_max})

  # return the reactiveValues
  return(filtering)
}
