#'
#' @import shinyWidgets
#' @import RColorBrewer
#' @import scales
#' @import esquisse
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
  list(`Brewer [divergent]`=list(`brewer:RdYlGn`=brewer_pal(palette='RdYlGn')(8),
                                 `brewer:PRGn`=brewer_pal(palette='PRGn')(8),
                                 `brewer:RdYlBu`=brewer_pal(palette='RdYlBu')(8),
                                 `brewer:RdGy`=brewer_pal(palette='RdGy')(8),
                                 `brewer:RdBu`=brewer_pal(palette='RdBu')(8),
                                 `brewer:PuOr`=brewer_pal(palette='PuOr')(8),
                                 `brewer:PiYG`=brewer_pal(palette='PiYG')(8),
                                 `brewer:BrBG`=brewer_pal(palette='BrBG')(8),
                                 `brewer:Spectral`=brewer_pal(palette='Spectral')(8)),
       `Brewer [sequential]`=list(`brewer:Blues`=brewer_pal(palette='Blues')(8),
                                  `brewer:BuPu`=brewer_pal(palette='BuPu')(8),
                                  `brewer:GnBu`=brewer_pal(palette='GnBu')(8),
                                  `brewer:Greens`=brewer_pal(palette='Greens')(8),
                                  `brewer:Greys`=brewer_pal(palette='Greys')(8),
                                  `brewer:Oranges`=brewer_pal(palette='Oranges')(8),
                                  `brewer:Purples`=brewer_pal(palette='Purples')(8),
                                  `brewer:RdPu`=brewer_pal(palette='RdPu')(8),
                                  `brewer:Reds`=brewer_pal(palette='Reds')(8),
                                  `brewer:YlGn`=brewer_pal(palette='YlGn')(8),
                                  `brewer:YlGnBu`=brewer_pal(palette='YlGnBu')(8)),
       `Viridis [sequential]`=list(`viridis:magma`=viridis_pal(option='magma')(8),
                                   `viridis:plasma`=viridis_pal(option='plasma')(8),
                                   `viridis:inferno`=viridis_pal(option='inferno')(8),
                                   `viridis:viridis`=viridis_pal(option='viridis')(8))) %>%
    palettePicker(inputId=NS(id, 'predefined_palette'), label='Select a palette', 
                  selected='brewer:YlGnBu', textColor=rgb(red=0, green=0, blue=0, alpha=0),
                  pickerOpts=list(`live-search`=FALSE, size=10)) -> spectrum_picker

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
