#'
#' @details \code{tags$style(type='text/css', '.irs-slider.from, .irs-slider.to {visibility: hidden !important;}')} to remove hide the sliders
components_picker_range.ui <- function(id, label='Principle components')
  sliderInput(inputId=NS(id, 'components_slider'), label=label, min=0, max=1, step=1, value=c(0,1), dragRange=TRUE)

#'
#' 
component_picker_slider.ui <- function(id, label='Principle component')
  sliderInput(inputId=NS(id, 'component_slider'), label=label, min=0, max=1, step=1, value=1, dragRange=FALSE)

#'
#' 
component_picker.ui <- function(id, label='Principle component')
  pickerInput(inputId=NS(id, 'component_picker'), label=label, choices=NULL, selected=NULL, multiple=FALSE)


#'
#' 
components_selector.server <- function(input, output, session, seurat, picked_reduction, range_size=5) {
  components <- reactiveValues(alpha=1)

  # react to a selection of component(s)
  ## from the range slider
  observe({req(input$components_slider); components$range <- input$components_slider})

  ## from the picker slider
  observe({req(input$component_slider); components$picked <- input$component_slider})

  ## from the picker
  observe({req(input$component_picker); components$picked <- input$component_picker})

  # react to a change in reduction method
  observe({
    req(seurat$object)
    req(picked_reduction$method)

    # get the number of components available in this reduction
    Stdev(object=seurat$object, reduction=picked_reduction$method) %>%
      seq_along() %>%
      purrr::set_names() -> principle_components
    min_value <- min(principle_components)
    max_value <- max(principle_components)

    # update ui elements
    ## update the range slider
    updateSliderInput(session=session, inputId='components_slider', min=-1)
    updateSliderInput(session=session, inputId='components_slider',
                      min=min_value, max=max_value, value=c(min_value,min_value+range_size-1))

    ## update the picker slider
    updateSliderInput(session=session, inputId='component_slider', min=-1)
    updateSliderInput(session=session, inputId='component_slider',
                      min=min_value, max=max_value, value=max_value)

    ## update the picker
    updatePickerInput(session=session, inputId='component_picker',
                      choices=principle_components, selected=min_value)})

  # return the reactiveValues list
  components
}
