seurat_object_options.ui <- function(id, seurat) {
  selectizeInput_defaults <- list(choices=NULL, selected=NULL, multiple=FALSE)
  list(inputId=NS(id, 'n_features_picker'), label='Features') %>% modifyList(x=selectizeInput_defaults) %>% do.call(what=selectizeInput) -> n_features_picker
  list(inputId=NS(id, 'n_umi_picker'), label='UMI') %>% modifyList(x=selectizeInput_defaults) %>% do.call(what=selectizeInput) -> n_umi_picker
  list(inputId=NS(id, 'proportion_mt_picker'), label='Mitochondrial proportion') %>% modifyList(x=selectizeInput_defaults) %>% do.call(what=selectizeInput) -> proportion_mt_picker

  # return ui element(s)
  tagList(tags$label('Configure metadata columns'),
          n_features_picker, n_umi_picker, proportion_mt_picker)
}

seurat_object_options.server <- function(input, output, session, seurat) {
}