#'
#' 
provenace_picker.ui <- function(id, seurat) {
  # get the analysis steps
  seurat$object@misc$provenance %>%
    names() %>%
    purrr::set_names() %>%
    purrr::set_names(function(x) str_c(seq_along(x), x, sep=': ')) -> analysis_steps

  selectInput(inputId=NS(id, 'analysis_step'), label='Analysis step',
              choices=analysis_steps, selected=analysis_steps[1])
}

#'
#' 
provenace_picker.server <- function(input, output, session, seurat) {
  provenace_picker <- reactiveValues(step=1, script='R script ...')

  observeEvent(eventExpr=input$analysis_step, label='provenace_picker/analysis_step', handlerExpr={
    provenace_picker$step <- input$analysis_step
    provenace_picker$script <- seurat$object@misc$provenance[[input$analysis_step]]})

  # return the reactiveValues list
  provenace_picker
}
