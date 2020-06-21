#'
#' 
point_size.ui <- function(id, label='Point size')
  sliderInput(inputId=NS(id, 'size_slider'), label=label, min=0.1, max=3, step=0.1, value=0.6)

#'
#'
point_size.server <- function(input, output, session, ...) {
  point_size <- reactiveValues(size=0.6)

  observe({req(input$size_slider); point_size$size <- input$size_slider})

  # return the reactiveValues list
  point_size
}
