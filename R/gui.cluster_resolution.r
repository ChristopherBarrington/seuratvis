#'
#' 
cluster_picker.ui <- function(id, seurat, resolution, picker, label_switch) {
  if(!missing(resolution) && is.logical(resolution) && resolution) resolution <- list()
  if(!missing(picker) && is.logical(picker) && picker) picker <- list()
  if(!missing(label_switch) && is.logical(label_switch) && label_switch) label_switch <- list()

  # get the cluster columns from the metadata
  tagList(if(is.list(resolution))
            list(inputId=NS(id, 'resolution_picker'), label='Cluster resolutions',
                 choices=seurat$cluster_resolutions, selected=preferred_choice(seurat$cluster_resolutions, 'seurat_clusters'), multiple=FALSE) %>%
              modifyList(val=resolution) %>%
              do.call(what=selectInput),
          if(is.list(picker))
            list(inputId=NS(id, 'ident_picker'), label='Cluster selection',
                 choices=seurat$all_idents$seurat_clusters, selected=seurat$all_idents$seurat_clusters, multiple=TRUE,
                 options=list(`actions-box`=TRUE, header='Cluster selection', title='Cluster selection',
                              `selected-text-format`='count>5', `count-selected-text`='{0} cluster(s)')) %>%
              modifyList(val=picker) %>%
              do.call(what=pickerInput),
          if(is.list(label_switch))
            list(inputId=NS(id, 'label_clusters'), label='Cluster labels',
                 value=TRUE, right=TRUE, status='success') %>%
              modifyList(val=label_switch) %>%
              do.call(what=materialSwitch))
}

#'
#' @import gtools
#' 
cluster_picker.server <- function(input, output, session, seurat, ...) {
  clusters <- reactiveValues()

  # react to a resolution being picked
  observeEvent(eventExpr=input$resolution_picker, label='cluster_picker/resolution_picker', handlerExpr={
    clusters$idents <- FetchData(object=seurat$object, vars=input$resolution_picker) %>% unlist()
    clusters$distinct_idents <- unique(clusters$idents) %>% mixedsort()
    updatePickerInput(session=session, inputId='ident_picker', choices=seurat$all_idents[[input$resolution_picker]], selected=seurat$all_idents[[input$resolution_picker]])})

  # react to idents being picked
  observeEvent(eventExpr=input$ident_picker, label='cluster_picker/ident_picker', handlerExpr={
    clusters$picked_idents <- input$ident_picker})

  # react to idents being picked
  observeEvent(eventExpr=input$label_clusters, label='cluster_picker/label_clusters', handlerExpr={
    clusters$label_clusters <- input$label_clusters})

  # return the reactiveValues list
  return(clusters)
}
