#' Displya R script used for analysis
#' 
#' Provides a selector and formatted code display of scripts
#' 
#' @param id unique name of the ekement
#' 
#' @imports shinyAce
#' 
#' @rdname provenance_step_viewer
#' 
provenance_step_viewer.ui <- function(id) {
  sprintf(fmt='### %s-provenance_step_viewer.ui', id) %>% message()

  module <- 'provenance_step_viewer'

  # make unique id for this object
  ns <- NS(namespace=id)
  module_ns <- ns(id=module)

  # create an environment in the seuratvis namespace
  e <- new.env()
  e$id <- id
  assign(x=module_ns, val=e, envir=module_environments)

  # track the re-used UI elements in each namespace
  get0(env=ui_element_ids.env, x=NS(namespace=module, id='analysis_step')) %>%
    append(ns(id='analysis_step')) %>%
    assign(env=ui_element_ids.env, x=NS(namespace=module, id='analysis_step'))

  # record the server(s) to call
  get0(env=module_servers_to_call, x=id) %>% append(sprintf(fmt='%s.server', module)) %>% assign(env=module_servers_to_call, x=id)

  # make ui elements
  ## analysis steps drop down box
  selectizeInput(inputId=ns(id='analysis_step'), label='Analysis step', choices=NULL, selected=NULL) -> analysis_step_picker

  ## output text box
  aceEditor(outputId='provenance_text', placeholder='R script',
            mode='r', theme='xcode',
            tabSize=2, useSoftTabs=TRUE, wordWrap=TRUE,
            showInvisibles=FALSE, highlightActiveLine=TRUE) -> text_output_box

  # return ui element(s)
  tagList(analysis_step_picker,
          text_output_box)
}

#' React to analysis step choice
#' 
#' @details
#' Picks analysis step from the \code{misc@provenance} slot of the Seurat object.
#' 
#' @imports shinyAce
#' @importFrom purrr set_names
#' 
#' @rdname provenance_step_viewer
#' 
provenance_step_viewer.server <- function(input, output, session) {
  session$ns('') %>% sprintf(fmt='### %sprovenance_step_viewer.server') %>% message()

  # get environments containing variables to run/configure this object
  collect_environments(id=parent.frame()$id, module='provenance_step_viewer') # provides `seuratvis_env`, `server_env` and `module_env`
  session_server <- get(x='session', env=server_env)

  # render the R code
  observeEvent(input$analysis_step, {
    # make sure these elements are defined
    req(seurat_object.reactions$seurat)

    # send a message
    sprintf(fmt='### %sprovenance_step_viewer.server-observeEvent-input$analysis_step [%s]', session$ns(''), input$analysis_step) %>% message()

    # create variables for shorthand
    r_script <- seurat_object.reactions$seurat@misc$provenance[[input$analysis_step]]
   
    # update the ui element(s)
    updateAceEditor(session=session_server, editorId='provenance_text', value=r_script)})

  # update UI when Seurat object is loaded
  observeEvent(eventExpr=seurat_object.reactions$seurat, handlerExpr={
    # send a message
    sprintf(fmt='### %sprovenance_step_viewer.server-observeEvent-seurat_object.reactions$seurat [%s]', session$ns(''), seurat_object.reactions$formatted.project.name) %>% message()

    # create variables for shorthand
    seurat <- seurat_object.reactions$seurat

    # get the analysis steps
    seurat@misc$provenance %>%
      names() %>%
      purrr::set_names() %>%
      purrr::set_names(function(x) str_c(seq_along(x), x, sep=': ')) -> analysis_steps
    picked_analysis_step <- analysis_steps[1]

    # update the ui element(s)
    ## metadata names dropdown box
    updateSelectizeInput(session=session, inputId='analysis_step',
                         choices=analysis_steps, selected=picked_analysis_step)})
}
