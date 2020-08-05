#'
#' 
add_feature_module_score.ui <- function(id) {
  ns <- NS(id)
  ui_width <- '400px'
  n_control_input <- numericInput(inputId=ns('n_control'), label='Number of control features', value=100, width=ui_width)
  module_name_input <- textInput(inputId=ns('module_name'), label='Name of module', placeholder='My module', width=ui_width, value='my_mod')
  textAreaInput(inputId=ns('feature_names'), label='Feature names', value='SOX2 ELAVL3',
    placeholder='Feature names separated by comma or whitespace (including newlines)',
    resize='vertical', rows=10, cols=40, width=ui_width) -> feature_names_input
  actionBttn(inputId=ns('add_score_button'), label='add module', style='stretch', color='warning') -> action_button

  tagList(module_name_input, n_control_input, feature_names_input, action_button)
}

#'
#' @import Seurat
#' 
add_feature_module_score.server <- function(input, output, session, seurat) {
  observeEvent(label='add_feature_module_score/add score', eventExpr=input$add_score_button, handlerExpr={
    req(seurat$gene_modules) #! TODO: these are features, not genes!
    req(seurat$gene_module_scores)
    req(input$module_name)
    req(input$feature_names)

print('add_feature_module_score/add score')
print(input$module_name)
print(input$feature_names)
    # if this is the first module and non were included, remove the dummy module
    if(seurat$gene_modules[1]=='dummy_module') {
      seurat$gene_modules <- c()
      seurat$gene_module_scores %<>% select(NULL) # does this work??!!
    }

    # get the list of feature names
    feature_names <- str_split(string=input$feature_names, pattern=',|\\s')
print(feature_names)
    # calculate the scores
    Seurat::AddModuleScore(object=seurat$object, features=feature_names, name='AddFeatureModuleServerGeneModule', ctrl=input$n_control) %>%
      slot(name='meta.data') %>%
      select(matches('^AddFeatureModuleServerGeneModule')) %T>% (function(x) print(head(x)))  %>%
      purrr::set_names(input$module_name) %T>% (function(x) print(head(x))) %>%
      cbind(seurat$gene_module_scores) -> seurat$gene_module_scores

    seurat$gene_modules %<>% c(input$module_name)

    # reset the ui elements
    reset(session$ns('module_name'))
    reset(session$ns('feature_names'))
  })
}







    # seurat$gene_modules <- c(dummy_module='normal distribution')
    # seurat$gene_module_scores <- seurat$metadata %>% select(NULL) %>% mutate(dummy_module=rnorm(n()))
    # gm_regex <- input$gene_modules_regex_text

    # if(!is.null(seurat$object@misc$gene_modules)) {
    #   seurat$gene_modules <- seurat$object@misc$gene_modules %>% purrr::set_names(str_remove, pattern=regex(gm_regex))
    #   seurat$gene_module_scores <- dplyr::select(seurat$metadata, matches(gm_regex)) %>% purrr::set_names(str_remove, pattern=regex(gm_regex))
    # }
    # seurat$gene_modules_regex <- gm_regex})



