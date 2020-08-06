#'
#' @import shinyWidgets
#' @import RColorBrewer
#'
colour_picker.ui <- function(id, label='Plot colours', include=c(Low='low', Mid='mid', High='high', `Plot background`='background'), gradient_selector_swictch=FALSE) {
  list(list('black', 'white', 'blanchedalmond', 'steelblue', 'forestgreen', 'linen', 'ivory', 'darkviolet', 'white'),
       as.list(brewer_pal(palette='Blues')(9)),
       as.list(brewer_pal(palette='Greens')(9)),
       as.list(brewer_pal(palette='Spectral')(11)),
       as.list(brewer_pal(palette='Dark2')(8))) -> possible_colours
    
  list(low='linen', mid='ivory', high='darkviolet', background='white') -> default_colours

  # make the pickers for each `include`
  ## define a set of options common to all pickers
  list(choices=possible_colours,
       flat=FALSE,
       options=list(`show-palette-only`=FALSE,
                    `toggle-palette-only`=FALSE,
                    `show-alpha`=FALSE,
                    `hide-after-palette-select`=FALSE,
                    `clickout-fires-change`=TRUE),
       update_on='move',
       width='50%') -> default_picker_options

  ## make a picker
  lapply(names(include), function(label)
    list(inputId=NS(id, include[[label]]),
         label=label,
         selected=default_colours[[include[[label]]]]) %>%
    modifyList(x=default_picker_options) %>%
    do.call(what=spectrumInput)) -> pickers

  # make a dropdown of pre-defined spectra
  list(`Brewer [d]`=c(`Spectral`='brewer:Spectral',
                      `RdYlGn`='brewer:RdYlGn',
                      `PRGn`='brewer:PRGn',
                      `RdYlBu`='brewer:RdYlBu',
                      `RdGy`='brewer:RdGy',
                      `RdBu`='brewer:RdBu',
                      `PuOr`='brewer:PuOr',
                      `PiYG`='brewer:PiYG',
                      `BrBG`='brewer:BrBG'),
       `Viridis [c]`=c(`magma`='viridis:magma',
                       `plasma`='viridis:plasma',
                       `inferno`='viridis:inferno',
                       `viridis`='viridis:viridis')) %>% print() %>%
    pickerInput(inputId=NS(id,'predefined_palette'), label='Select a palette',
                multiple=FALSE, options=list(`live-search`=TRUE, title='Colour spectra')) -> spectrum_picker

  # rerurn the ui
  tagList(tags$style(type='text/css', '.sp-dd {display:none !important} .sp-preview {width:100% !important}'),
          tags$label(label),
          pickers,
          spectrum_picker)
}

#'
#'
colour_picker.server <- function(input, output, session) {
  colours <- reactiveValues(low='linen', mid='ivory', high='darkviolet', background='white')

  # react to the colours being picked
  ## low
  observe({req(input$low); colours$low <- input$low})

  ## mid
  observe({req(input$mid); colours$mid <- input$mid})

  ## high
  observe({req(input$high); colours$high <- input$high})

  ## spectrum
  observe({req(input$predefined_palette); colours$palette <- input$predefined_palette %>% str_split(pattern=':') %>% pluck(1)})

  ## plot background
  observe({req(input$background); colours$background <- input$background})

  # return the reactives list
  colours
}
