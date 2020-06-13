#' Select a principle component from a dataset PCA reduction
#' 
#' Provides a method to select principle component(s) from a PCA reduction
#' 
#' @param id unique name of the element
#' @param label text label of the element
#' 
#' @examples
#' 
#' @rdname principle_component_picker
#' 
principle_component_picker.ui <- function(id, label='Principle component selection', include_picker=TRUE, include_slider=TRUE) {
  sprintf(fmt='### %s-principle_component_picker.ui', id) %>% message()

  module <- 'principle_component_picker'

  # make unique id for this object
  ns <- NS(namespace=id)
  module_ns <- ns(id=module)

  # create an environment in the seuratvis namespace
  e <- new.env()
  e$id <- id
  assign(x=module_ns, val=e, envir=module_environments)

  # track the re-used UI elements in each namespace
  get0(env=ui_element_ids.env, x=NS(namespace=module, id='principle_component_picker')) %>%
    append(ns(id='principle_component_picker')) %>%
    assign(env=ui_element_ids.env, x=NS(namespace=module, id='principle_component_picker'))

  # record the server(s) to call
  get0(env=module_servers_to_call, x=id) %>% append(sprintf(fmt='%s.server', module)) %>% assign(env=module_servers_to_call, x=id)

  # make ui elements
  ## principle components drop down box
  selectizeInput(inputId=ns(id='principle_component_picker'), label='Principle component', choices=NULL, selected=NULL) -> principle_component_picker

  ## slider to limit range of components
  sliderInput(inputId=ns(id='principle_components_slider'), label='Principle component',
              min=0, max=1, step=1, value=-Inf, dragRange=FALSE) -> principle_components_slider

  # return ui element(s)
  # tagList(dropdown, metadata_switch, value_range)
  tagList(if(include_picker) principle_component_picker,
          if(include_slider) principle_components_slider)
}

#' React to a principle component choice
#' 
#' @details
#' Updates the Seurat object reactive values.
#' 
#' @importFrom purrr set_names
#' 
#' @rdname principle_component_picker
#' 
principle_component_picker.server <- function(input, output, session) {
  session$ns('') %>% sprintf(fmt='### %sprinciple_component_picker.server') %>% message()

  # get environments containing variables to run/configure this object
  collect_environments(id=parent.frame()$id, module='principle_component_picker') # provides `seuratvis_env`, `server_env` and `module_env`
  session_server <- get(x='session', env=server_env)

  # react to the component selection
  ## if a component is selected, copy it to the reactive
  observeEvent(eventExpr={input$principle_component_picker}, ignoreInit=TRUE, handlerExpr={
    # make sure these elements are defined
    req(selections.rv[[session$ns('picked_component')]])

    # send a message
    sprintf(fmt='### %sprinciple_component_picker.server-observeEvent-input$principle_component_picker[%s]', session$ns(''), input$principle_component_picker) %>% message()
    
    selections.rv[[session$ns('picked_component')]] <- input$principle_component_picker})

  # react to the colour scale limits
  observeEvent(eventExpr=input$principle_components_slider, handlerExpr={
    # make sure these elements are defined
    req(selections.rv[[session$ns('value_range_limits')]])

    # send a message
    sprintf(fmt='### %sprinciple_component_picker.server-observeEvent-input$principle_components_slider [%s]', session$ns(''), str_c(input$principle_components_slider, collapse=',')) %>% message()

    # update the reactive
    selections.rv[[session$ns('picked_components_range')]] <- input$principle_components_slider})

  # update UI when reduction type is changed
  observeEvent(eventExpr=input$reduction_method_picker, handlerExpr={
    # make sure these elements are defined
    req(seurat_object.reactions$seurat)
    req(input$reduction_method_picker)

    # send a message
    sprintf(fmt='### %sprinciple_component_picker.server-observeEvent-selections.rv[[session$ns(selected_reduction_method)]] [%s]', session$ns(''), seurat_object.reactions$formatted.project.name) %>% message()

    # create variables for shorthand
    seurat <- seurat_object.reactions$seurat
    reduction_name <- input$reduction_method_picker

    # get the number of components available in this reduction
    Stdev(object=seurat, reduction=reduction_name) %>%
      seq_along() %>%
      purrr::set_names() -> principle_component_picker_options
    min_value <- min(principle_component_picker_options)
    max_value <- max(principle_component_picker_options)

    # update the ui element(s)
    ## principle components drop down box
    updateSelectizeInput(session=session, inputId='principle_component_picker',
                         choices=principle_component_picker_options, selected=min_value)
    
    ## slider to limit range of components
    updateSliderInput(session=session, inputId='principle_components_slider',
                      min=min_value, max=max_value, value=min_value+5)

    # update the reactive
    selections.rv[[session$ns('picked_component')]] <- min_value
    selections.rv[[session$ns('picked_components_range')]] <- c(min_value, max_value)})
}
