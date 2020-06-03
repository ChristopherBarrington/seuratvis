#' Specify maximum mitochondrial proportion
#' 
#' Provides a way to define limit on proportion of mitochondrial UMI per cell
#' 
#' @param id unique name of the element
#' @param label text label of the element
#' 
#' @examples
#' 
#' \dontrun{
#' ui <- fluidPage(percent_mt_per_cell_filter.ui(id='page_name'))
#' server <- function(input, output, session) {
#'   callModule(percent_mt_per_cell_filter.server, id='page_name')
#' }}
#' 
#' @rdname percent_mt_per_cell_filter
#' 
percent_mt_per_cell_filter.ui <- function(id, label='Proportion Mt') {
  sprintf(fmt='### %s-percent_mt_per_cell_filter.ui', id) %>% message()

  module <- 'percent_mt_per_cell_filter'

  # make unique id for this object
  ns <- NS(namespace=id)
  module_ns <- ns(id=module)

  # create an environment in the seuratvis namespace
  e <- new.env()
  e$id <- id
  assign(x=module_ns, val=e, envir=module_environments)


  # record the server(s) to call
  module_environments$percent_mt_per_cell_filters$id %<>% c(id) # keep track so other modules can update
  get0(env=module_servers_to_call, x=id) %>% append(sprintf(fmt='%s.server', module)) %>% assign(env=module_servers_to_call, x=id)

  # make ui elements
  div(tags$h6('Maximum', style='display: inline;'),
    numericInput(inputId=ns(id='max_percent_mt'), label=NULL, value=0, step=0.1, width='100%')) -> high_ui

  # return ui element(s)
  tagList(tags$label(label), br(), high_ui)
}

#' React to maximum proportion mitochondria per cell threshold
#' 
#' @details
#' Updates the Seurat object reactive values.
#' 
#' @rdname percent_mt_per_cell_filter
#'
percent_mt_per_cell_filter.server <- function(input, output, session) {
  session$ns('') %>% sprintf(fmt='### %spercent_mt_per_cell_filter.server') %>% message()

  # get environments containing variables to run/configure this object
  collect_environments(id=parent.frame()$id, module='percent_mt_per_cell_filter') # provides `seuratvis_env`, `server_env` and `module_env`

  # react to the maximum input element
  observeEvent(eventExpr=input$max_percent_mt, handlerExpr={
    message('### percent_mt_per_cell_filter.server-observeEvent-input$max_percent_mt')

    # update the reactive
    value <- input$max_percent_mt

    if(value != sprintf(fmt='%.1f', value)) # if the value is not already 1dp formatted, reformat it
      value %<>% add(0.05) %>% round(digits=1)

    filtering_parameters.reactions$max_percent_mitochondria <- value})

  # react to the initialisation of the reference value
  observeEvent(eventExpr=seurat_object.reactions$proportion_mt_values_max, handlerExpr={
    message('### percent_mt_per_cell_filter.server-observeEvent-seurat_object.reactions$proportion_mt_values_max')

    # create variables for shorthand
    high <- seurat_object.reactions$proportion_mt_values_max %>% add(0.05) %>% round(digits=1)

    # update the ui element(s)
    updateNumericInput(session=session, inputId='max_percent_mt', value=high, min=0, max=high)})
}
