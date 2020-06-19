#' Display a colour selection panel
#' 
#' Create UI elements for colours and a palette type switch
#' 
#' @param id unique name of the element
#' @param selectors list of lists with arguments to \code{colourInput}; requires \code{inputId}, \code{label} and best include \code{value} too
#' @param include_full boolean whether to create a square/limited palette switch
#' 
#' @details
#' For every selector in the list, a new element is created in the same order as provided. If a palette type switch is required, it is created after the selectors.
#' 
#' @examples
#' 
#' \dontrun{
#' ui <- fluidPage(number_of_clusters_text_box.ui(id='page_name'))
#' server <- function(input, output, session) {
#'   callModule(module=update_palette_type.server, id='page_name')
#'   callModule(module=add_to_colour_palette.server, id='page_name')
#' }}
#' 
#' @import colourpicker
#' 
#' @rdname colour_palette
#' 
colour_palette.ui <- function(id, label='Feature value colours', selectors=list(), include_full=FALSE) {
  sprintf(fmt='### %s-colour_palette.ui', id) %>% message()

  module <- 'colour_palette'

  # make unique id for this object
  module_ns <- ns <- NS(namespace=id, id=module)
  ns <- NS(namespace=id)

  # record the server(s) to call
  get0(env=module_servers_to_call, x=id) %>% append(sprintf(fmt='%s.server', c('update_palette_type', 'add_to_colour_palette'))) %>% assign(env=module_servers_to_call, x=id)

  # colour label names should be one of 'low', 'mid' or 'high'
  if(any(! sapply(selectors, pluck, 'inputId') %in% c('low','mid','high')))
    stop("!!! inputID to colour_palette.ui show be one of 'low', 'mid' or 'high'")

  # for each selector, make a list UI element arguments
  defaults <- list(value=sample(x=default_colour_palette(), size=1), showColour='both', palette='limited', allowedCols=default_colour_palette(), allowTransparent=FALSE, returnName=TRUE)
  selectors %<>%
    lapply(function(x) {x$inputId %<>% ns(); x})

  # get inputIDs
  selector_inputIds <- sapply(selectors, function(x) x$inputId)

  # make colour selector UIs
  selectors %<>%
    lapply(modifyList, x=defaults) %>%
    lapply(do.call, what=colourInput)

  # if full palette is chosen, add the switch
  full_palette_switch <- materialSwitch(inputId=ns(id='full'), label='Show full palette?', value=FALSE, right=TRUE, status='success')
  if(include_full)
    selectors %<>% append(list(full_palette_switch))

  # create an environment in the seuratvis namespace
  e <- new.env()
  e$id <- id
  e$selector_inputIds <- selector_inputIds
  e$full_palette_inputId <- ns(id='full')
  e$include_full_palette_switch <- include_full
  assign(x=module_ns, val=e, envir=module_environments)

  # return ui element(s)
  tagList(tags$label(label), br(), selectors)
}

#' Update the colour selectors palette style
#' 
#' @rdname colour_palette
#' 
update_palette_type.server <- function(input, output, session, ...) {
  session$ns('') %>% sprintf(fmt='### %supdate_palette_type.server') %>% message()

  # get environments containing variables to run/configure this object
  collect_environments(id=parent.frame()$id, module='colour_palette') # provides `seuratvis_env`, `server_env` and `module_env`
  input <- get(x='input', env=server_env)
  session <- get(x='session', env=server_env)

  # react to the palette switch
  observe({
    if(!module_env$include_full_palette_switch) # limited palette doesn't work without this ... ?
      return(NULL)

    palette_type <- ifelse(input[[module_env$full_palette_inputId]], 'square', 'limited')

    for(inputId in module_env$selector_inputIds)
      updateColourInput(session=session, inputId=inputId, palette=palette_type, value=input[[inputId]], allowedCols=colour_palette())})
}

#' Add a new colour to the colour palette
#' 
#' @rdname colour_palette
#' 
add_to_colour_palette.server <- function(input, output, session, ...) {
  session$ns('') %>% sprintf(fmt='### %sadd_to_colour_palette.server') %>% message()

  # get environments containing variables to run/configure this object
  collect_environments(id=parent.frame()$id, module='colour_palette') # provides `seuratvis_env`, `server_env` and `module_env`
  input_server <- get(x='input', env=server_env)

  # react to the chosen colour by adding it to the colour palette reactive values
  observe({
    for(inputId in module_env$selector_inputIds) {
      x <- input_server[[inputId]]
      if(x!='' && !startsWith(x=x, prefix='#')) # if x is not in hex format
        x %<>% gplots::col2hex()
      unique(c(colour_palette(), x)) %>% colour_palette()
    }})
}

#' Define default colour palette
#' 
default_colour_palette <- function()
  c('#000000', '#4D4D4D', '#666666', '#7F7F7F', '#999999',
    '#B3B3B3', '#E5E5E5', '#FFFFFF', '#FAF0E6', '#27408B',
    '#000080', '#0000FF', '#1E90FF', '#63B8FF', '#97FFFF',
    '#00FFFF', '#00868B', '#008B45', '#458B00', '#008B00',
    '#00FF00', '#7FFF00', '#54FF9F', '#00FF7F', '#7FFFD4',
    '#8B4500', '#8B0000', '#FF0000', '#FF6A6A', '#FF7F00',
    '#FFFF00', '#FFF68F', '#F4A460', '#551A8B', '#8B008B',
    '#8B0A50', '#9400D3', '#FF00FF', '#FF1493', '#E066FF')

#' A reactive to store the colours
#' 
#' This is probably not a good place to keep this
#' 
colour_palette <- reactiveVal(default_colour_palette()) # should this be here?


