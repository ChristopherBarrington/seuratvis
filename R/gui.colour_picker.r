#'
#' @import shinyWidgets
#' @import RColorBrewer
#'
colour_picker.ui <- function(id, label='Plot colours', include=c(Low='low', Mid='mid', High='high', `Plot background`='background')) {

  # define a list of colours that can be picked
  list(list('black', 'white', 'ivory', 'linen', 'blanchedalmond'),
      as.list(brewer_pal(palette='Blues')(9)),
      as.list(brewer_pal(palette='Greens')(9)),
      as.list(brewer_pal(palette='Greys')(9)),
      as.list(brewer_pal(palette='Oranges')(9)),
      as.list(brewer_pal(palette='Purples')(9)),
      as.list(brewer_pal(palette='Reds')(9)),
      as.list(brewer_pal(palette='PiYG')(11)),
      as.list(brewer_pal(palette='RdBu')(11)),
      as.list(brewer_pal(palette='Spectral')(11))) -> possible_colours
   
  list(low='linen', mid='ivory', high='#5E4FA2', background='white') -> default_colours

  # make the pickers for each `include`
  ## define a set of options common to all pickers
  list(choices=possible_colours,
      flat=FALSE,
      options=list(`show-palette-only`=TRUE,
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
  list(`Brewer [divergent]`=c(`Red-Yellow-Green`='brewer:RdYlGn',
                              `Purple-Green`='brewer:PRGn',
                              `Red-Yellow-Blue`='brewer:RdYlBu',
                              `Red-Grey`='brewer:RdGy',
                              `Red-Blue`='brewer:RdBu',
                              `Orange-Purple`='brewer:PuOr',
                              `Pink-Green`='brewer:PiYG',
                              `Brown-Teal`='brewer:BrBG',
                              `Spectral`='brewer:Spectral'),
       `Brewer [sequential]`=c(`Blues`='brewer:Blues',
                               `Blue-Purple`='brewer:BuPu',
                               `Green-Blue`='brewer:GnBu',
                               `Greens`='brewer:Greens',
                               `Greys`='brewer:Greys',
                               `Oranges`='brewer:Oranges',
                               `Purples`='brewer:Purples',
                               `Red-Purple`='brewer:RdPu',
                               `Reds`='brewer:Reds',
                               `Yellow-Green`='brewer:YlGn',
                               `Yellow-Green-Blue`='brewer:YlGnBu'),
       `Viridis [sequential]`=c(`magma`='viridis:magma',
                                `plasma`='viridis:plasma',
                                `inferno`='viridis:inferno',
                                `viridis`='viridis:viridis')) %>%
    pickerInput(inputId=NS(id, 'predefined_palette'), label='Select a palette', selected='viridis:plasma',
                multiple=FALSE, options=list(`live-search`=FALSE, size=10), inline=FALSE) -> spectrum_picker

  # make a selection for gradient direction
  prettyToggle(inputId=NS(id, 'palette_direction'),
               label_on='Forward palette', icon_on=icon('sort-amount-up'), status_on='success', 
               label_off='Reverse palette', icon_off=icon('sort-amount-down'), status_off='danger',
               value=TRUE, outline=TRUE, plain=TRUE, animation='jelly') -> palette_direction

  # rerurn the ui
  tagList(tags$style(type='text/css', '.sp-dd {display:none !important} .sp-preview {width:100% !important}'),
          tags$label(label),
          pickers,
          spectrum_picker,
          palette_direction)
}

#'
#'
colour_picker.server <- function(input, output, session) {
  colours <- reactiveValues(low='linen', mid='ivory', high='#5E4FA2', background='white', palette=list('viridis','plasma'), direction=1)

  # react to the colours being picked
  ## low
  observe({req(input$low); colours$low <- input$low})

  ## mid
  observe({req(input$mid); colours$mid <- input$mid})

  ## high
  observe({req(input$high); colours$high <- input$high})

  ## palette
  observe({req(input$predefined_palette); colours$palette <- input$predefined_palette %>% str_split(pattern=':') %>% pluck(1)})

  ## palette direction
  observe({if(is.null(input$palette_direction)) return(NULL); colours$direction <- input$palette_direction %>% ifelse(1, -1)})

  ## plot background
  observe({req(input$background); colours$background <- input$background})

  # return the reactives list
  colours
}
