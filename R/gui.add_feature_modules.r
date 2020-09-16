#'
#' 
add_feature_module_score.ui <- function(id) {
  ns <- NS(id)
  ui_width <- '400px'

  n_control_input <- numericInput(inputId=ns('n_control'), label='Number of control features', value=100, width=ui_width) %>% hidden()
  module_name_input <- textInput(inputId=ns('module_name'), label='Name of module', placeholder='My module', width=ui_width)
  textAreaInput(inputId=ns('features_in_module'), label='Feature names',
    placeholder='Feature names separated by comma or whitespace (including newlines)',
    resize='vertical', rows=10, cols=40, width=ui_width) -> feature_names_input
  actionButton(inputId=ns('add_score_button'), label='Add the module', style='simple', color='default') -> action_button #! TODO: this should be an actionBttn

  tagList(module_name_input,
          n_control_input,
          feature_names_input,
          action_button)
}

#'
#' @import Seurat
#' 
add_feature_module_score.server <- function(input, output, session, seurat) {
  observeEvent(label='add_feature_module_score/add score', eventExpr=input$add_score_button, handlerExpr={
    req(seurat$gene_modules) #! TODO: these are features, not genes!
    req(seurat$gene_module_scores)
    req(input$module_name)
    req(input$features_in_module)

    # if this is the first module and non were included, remove the dummy module
    if(names(seurat$gene_modules)[1]=='dummy_module') {
      seurat$gene_modules <- list()
      seurat$gene_module_scores %<>% dplyr::select(NULL)
    }

    # if module name already exists, append to make it unique
    if(is.element(el=input$module_name, names(seurat$gene_modules))) {
      error_alert(session=session, title='Choose unique module name', text=sprintf(fmt='%s exists already!', input$module_name))
      return(NULL)
    }

    # get the list of feature names
    features_in_module <- str_split(string=input$features_in_module, pattern=',|\\s') %>% set_names(input$module_name)

    # calculate the scores
    Seurat::AddModuleScore(object=seurat$object, features=features_in_module, name='AddFeatureModuleServerGeneModule', ctrl=input$n_control) %>%
      slot(name='meta.data') %>%
      dplyr::select(matches('^AddFeatureModuleServerGeneModule')) %>%
      purrr::set_names(input$module_name) %>%
      cbind(seurat$gene_module_scores) -> seurat$gene_module_scores

    seurat$gene_modules %<>% append(features_in_module)

    # reset the ui elements
    reset(id='module_name')
    reset(id='features_in_module')
  })
}
