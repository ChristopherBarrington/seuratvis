#'
#' 
components_selector.ui <- function(id, label='Principle components')
  sliderInput(inputId=NS(id, 'components_slider'), label=label, min=0, max=1, step=1, value=c(40,50), dragRange=TRUE)

#'
#' 
component_picker.ui <- function(id, label='Principle component')
  pickerInput(inputId=NS(id, 'component_picker'), label=label, choices=NULL, selected=NULL, multiple=FALSE)


#'
#' 
components_selector.server <- function(input, output, session, seurat, picked_reduction, range_size=15) {
  components <- reactiveValues(alpha=1)

  # react to a selection of component(s)
  ## from the slider
  observe({req(input$components_slider); components$range <- input$components_slider})

  # from the picker
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
    ## update the slider
    updateSliderInput(session=session, inputId='components_slider', min=-1)
    updateSliderInput(session=session, inputId='components_slider',
                      min=min_value, max=max_value, value=c(min_value,min_value+range_size-1))

    ## update the picker
    updatePickerInput(session=session, inputId='component_picker',
                      choices=principle_components, selected=min_value)})

  # return the reactiveValues list
  components
}
