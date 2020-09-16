#'
#' 
dimension_reduction.ui <- function(id, seurat, regex='.*', label='Reduction method') {
  Reductions(seurat$object) %>%
    str_subset(pattern=regex(pattern=regex)) %>%
    str_subset(pattern=regex(pattern='3D', ignore_case=TRUE), negate=TRUE) -> reductions

  # do nothing if there are no reductions
  if(length(reductions)==0)
    return(NULL)

  selectInput(inputId=NS(id, 'picker'), label=label,
              choices=reductions, multiple=FALSE,
              selected=preferred_choice(x=reductions, preferences=c('umap','tsne','pca')))
}

#'
#' 
dimension_reduction.server <- function(input, output, session, seurat, regex='.*', ...) {
  reduction <- reactiveValues()

  # react to the input
  # observeEvent(eventExpr=input$picker, label='dimension_reduction/picker', handlerExpr={
  observe(label='dimension_reduction/picker', x={
    req(seurat$object)
    req(input$picker)

    # get the input
    picked <- input$picker

    # get the embeddings
    embeddings <- tryCatch(Embeddings(object=seurat$object, reduction=picked), error=function(...) return(data.frame(DIMRED_1=numeric(), DIMRED_2=numeric())))
    embeddings[,1:2] %>%
      as.data.frame() %>%
      set_names(c('DIMRED_1','DIMRED_2')) -> embeddings

    # update the reactiveValues
    reduction$embeddings <- embeddings
    reduction$method <- picked})

  # return the reactiveValues list
  return(reduction)
}
