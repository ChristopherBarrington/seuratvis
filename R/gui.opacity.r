#'
#' 
opacity.ui <- function(id, label='Opacity')
  sliderInput(inputId=NS(id, 'alpha_slider'), label=label, min=0.1, max=1, step=0.05, value=1)

#'
#' 
opacity.server <- function(input, output, session, ...) {
  opacity <- reactiveValues(alpha=1)

  observe({req(input$alpha_slider); opacity$alpha <- input$alpha_slider})

  # return the reactiveValues list
  opacity
}
