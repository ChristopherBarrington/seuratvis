#' 
#' 
provenance_step_viewer.ui <- function(id)
  aceEditor(outputId=NS(id, 'ace_editor'), placeholder='R script',
            mode='r', theme='xcode',
            tabSize=2, useSoftTabs=TRUE, wordWrap=TRUE,
            showInvisibles=FALSE, highlightActiveLine=TRUE)

#' 
#' 
provenance_step_viewer.server <- function(input, output, session, picked_provenance) {
  observeEvent(eventExpr=picked_provenance$script, label='provenance_step_viewer/script', handlerExpr={
    updateAceEditor(session=session, editorId='ace_editor', value=picked_provenance$script)
  })
}


ace_editor.ui <- function(id, ...)
  aceEditor(outputId=NS(id, 'ace_editor'), placeholder='R script',
            mode='r', tabSize=2, useSoftTabs=TRUE, wordWrap=TRUE,
            showInvisibles=FALSE, highlightActiveLine=FALSE, ...)

ace_editor.server <- function(input, output, session, display_text)
  observe(label='ace_editor/observe', x={
    updateAceEditor(session=session, editorId='ace_editor', value=display_text())})
