#' Make a PCA heatmap of a reduction
#' 
#' Creates a heatmap of features contributing to components
#' 
#' @param id unique name of the element
#' 
#' @rdname pca_heatmap
#' 
pca_heatmap.ui <- function(id) {
  sprintf(fmt='### %s-pca_heatmap.ui', id) %>% message()

  module <- 'pca_heatmap'

  # make unique id for this object
  ns <- NS(namespace=id)
  module_ns <- ns(id=module)

  # create an environment in the seuratvis namespace
  e <- new.env()
  e$id <- id
  assign(x=module_ns, val=e, envir=module_environments)

  # record the server(s) to call
  get0(env=module_servers_to_call, x=id) %>% append(sprintf(fmt='%s.server', module)) %>% assign(env=module_servers_to_call, x=id)

  # return ui element(s)
  tagList(plotOutput(outputId=ns(id='pca_heatmap')) %>% withSpinner())
}

#' Produce the ggplot object for a PCA heatmap
#' 
#' @rdname pca_heatmap
#'
pca_heatmap.server <- function(input, output, session, ...) {
  session$ns('') %>% sprintf(fmt='### %spca_heatmap.server') %>% message()

  # get environments containing variables to run/configure this object
  collect_environments(id=parent.frame()$id, module='pca_heatmap') # provides `seuratvis_env`, `server_env` and `module_env`
  id <- parent.frame()$id
  ns <- NS(namespace='pca_heatmap')

  # render the elbow plot
  renderPlot(expr={
    # make sure these elements are defined
    req(seurat_object.reactions$seurat)
    req(input$reduction_method_picker)
    req(input$principle_component_picker)

    # send a message
    session$ns('') %>% sprintf(fmt='### %spca_heatmap.server-renderPlot') %>% message('')
   
    # make variables for shorthand    
    seurat <- seurat_object.reactions$seurat
    reduction_name <- input$reduction_method_picker
    selected_component <- input$principle_component_picker

    # render the heatmap
    #! TODO: make this a nice ggplot
    if(!DefaultAssay(object=seurat[['pca']]) %in% Assays(seurat)) {
      ggplot()+aes()+annotate(geom='text', label='Nothing to see here', x=0, y=0)+theme_void()
    } else {
      DefaultAssay(object=seurat) <- DefaultAssay(object=seurat[[reduction_name]])
      DimHeatmap(object=seurat, reduction=reduction_name,
                 dims=as.numeric(selected_component),
                 disp.min=-2.5, disp.max=2.5,
                 cells=2000, balanced=TRUE, fast=FALSE)}}) -> output$pca_heatmap
}
