#'
#'
components_picker_range.ui <- function(id, label='Principle components') {
  # \code{tags$style(type='text/css', '.irs-slider.from, .irs-slider.to {visibility: hidden !important;}')} to remove hide the sliders
  sliderInput(inputId=NS(id, 'components_slider'), label=label, min=1, max=10, step=1, value=c(1,2), dragRange=TRUE) -> selector
  numericInput(inputId=NS(id, 'refresh_components_slider'), label='refresh_selectors', value=rnorm(1)) %>% hidden() -> refresher
  tagList(selector, refresher)
}

#'
#' 
component_picker_slider.ui <- function(id, label='Principle component') {
  sliderInput(inputId=NS(id, 'component_slider'), label=label, min=1, max=10, step=1, value=1, dragRange=FALSE) -> selector
  numericInput(inputId=NS(id, 'refresh_component_slider'), label='refresh_selectors', value=rnorm(1)) %>% hidden() -> refresher
  tagList(selector, refresher)
}

#'
#' 
component_picker.ui <- function(id, label='Principle component') {
  pickerInput(inputId=NS(id, 'component_picker'), label=label, choices=NULL, selected=NULL, multiple=FALSE) -> selector
  numericInput(inputId=NS(id, 'refresh_component_picker'), label='refresh_selectors', value=rnorm(1)) %>% hidden() -> refresher
  tagList(selector, refresher)
}

#'
#' 
components_selector.server <- function(input, output, session, seurat, picked_reduction, range_size=5) {
  components <- reactiveValues()

  # react to a selection of component(s)
  ## from the range slider
  observe({req(input$components_slider); components$range <- input$components_slider})

  ## from the picker slider
  observe({req(input$component_slider); components$picked <- input$component_slider})

  ## from the picker
  observe({req(input$component_picker); components$picked <- input$component_picker})

  # react to a change in reduction method
  observeEvent(eventExpr=c(input$refresh_components_slider, input$refresh_component_slider, input$refresh_component_picker), handlerExpr={
    req(seurat$n_principle_components)
    req(picked_reduction$method)

    # get the number of components available in this reduction
    min_value <- 1
    max_value <- seurat$n_principle_components[[picked_reduction$method]]

    # update ui elements
    ## update the range slider
    updateSliderInput(session=session, inputId='components_slider', max=max_value, value=c(1,2))

    ## update the picker slider
    updateSliderInput(session=session, inputId='component_slider', max=max_value, value=1)

    ## update the picker
    updatePickerInput(session=session, inputId='component_picker', choices=seq(max_value), selected=min_value)})

  # return the reactiveValues list
  components
}
