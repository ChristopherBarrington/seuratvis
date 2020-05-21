#' Display a table of available Seurat objects
#' 
#' Shows information on Seurats in accessible environments for selection
#' 
#' @param id unique name of the element
#' 
#' @details
#' A \code{datatable} is rendered with additional Seurat metadata. The selected row is loaded into the viewer.
#' 
#' @examples
#' 
#' \dontrun{
#' ui <- fluidPage(available_seurats.ui(id='page_name'))
#' server <- function(input, output, session) {
#'   callModule(module=available_seurats.server, id='page_name')
#'   callModule(module=load_a_seurat.server, id='page_name')
#' }}
#' 
#' @rdname available_seurats
#' 
available_seurats.ui <- function(id) {
  message('### available_seurats.ui')

  module <- 'available_seurats'

  # make unique id for this object
  ns <- NS(namespace=id, id=module)

  # create an environment in the seuratvis namespace
  e <- new.env()
  e$id <- id
  assign(x=ns, val=e, envir=module_environments)

  # return ui element(s)
  DT::dataTableOutput(outputId=ns)
}

#' Render the table of Seurat object and their metadata
#' 
#' @import DT
#' 
#' @rdname available_seurats
#' 
available_seurats.server <- function(input, output, session) {
  message('### available_seurats.server')

  # get environments containing variables to run/configure this object
  collect_environments(id=parent.frame()$id, module='available_seurats') # provides `seuratvis_env`, `server_env` and `module_env`

  # define output column names and order
  c(Project='project',
    Environment='environment',
    `Object name`='object',
    `Number of cells`='ncells',
    `Features in active assay`='nfeatures',
    `Filtered?`='filtered',
    `Integrated?`='integrated') -> column_order

  # make the `data.frame` of Seurat information
  server_env$available_seurat_objects %>%
    dplyr::select(-choiceName) %>%
    plyr::adply(.margins=1, function(params) {
      x <- eval(parse(text=params$choiceValue))
      data.frame(project=reformat_project_name(Project(x)),
                 environment=if_else(params$env=='globalenv()', 'RGlobal', as.character(params$env)),
                 ncells=ncol(x),
                 nfeatures=nrow(x),
                 dimensions=ifelse(is.null(x@misc$n_dimensions), -1, x@misc$n_dimensions),
                 filtered=!is.null(x@misc$cells_filtered) && x@misc$cells_filtered,
                 integrated=!is.null(x@misc$integrated_dataset) && x@misc$integrated_dataset)}) %>%
    mutate(dimensions=as.integer(dimensions),
           object=value) %>%
    select_at(vars(all_of(column_order), everything())) -> data_to_show

  # identify columns to format in the `datatable` `columnDefs` argument
  formatted_colnames <- c(names(column_order), colnames(data_to_show)[! colnames(data_to_show) %in% column_order])
  list(hide={seq(length(column_order)+1, ncol(data_to_show))},
       center={sapply(data_to_show, class) %>% str_which('logical')}) %>%
    lapply(subtract, e2=1) -> columnDef_targets

  # make and format the `datatable`
  DT::renderDataTable({
    data_to_show %>%
      DT::datatable(colnames=formatted_colnames,
                    rownames=FALSE,
                    options=list(columnDefs=list(list(visible=FALSE, targets=columnDef_targets$hide),
                                                 list(className='dt-center', targets=columnDef_targets$center)),
                                 ordering=FALSE,
                                 dom='ti',
                                 language=list(info='Found _TOTAL_ Seurat object(s)',
                                               infoEmpty='No Seurat objects found!')),
                    style='bootstrap4',
                    class='stripe',
                    selection=list(mode='single',
                                   selected=input$available_seurats_rows_selected)) %>%
      formatStyle(columns='ncells',
                  background=styleColorBar(data=c(0,max(data_to_show$ncells)), color='#3CB96A'),
                  backgroundSize='98% 50%',
                  backgroundRepeat='no-repeat',
                  backgroundPosition='center') %>%
      formatStyle(columns='nfeatures',
                  background=styleColorBar(data=c(0,max(data_to_show$nfeatures)), color='#3CB96A'),
                  backgroundSize='98% 50%',
                  backgroundRepeat='no-repeat',
                  backgroundPosition='center') %>%
      formatRound(columns=c('ncells', 'nfeatures'),
                  digits=0) %>%
      formatStyle(columns={sapply(data_to_show, class) %>% str_which('logical')},
                  color=styleEqual(levels=c(0,1), values=c('#B96A3C','#B93C8B'), default='orange'),
                  fontFamily='monospace',
                  fontWeight='bold')}) -> output$available_seurats
}

#' Load a seurat object into a reactive
#' 
#' @rdname available_seurats
#' 
load_a_seurat.server <- function(input, output, session) {
  message('### load_a_seurat.server')

  # get environments containing variables to run/configure this object
  collect_environments(id=parent.frame()$id, module='available_seurats') # provides `seuratvis_env`, `server_env` and `module_env`
  session <- get(x='session', env=server_env)

  # react to a row being selected in the available Seurats table
  observeEvent(eventExpr=input$available_seurats_rows_selected, handlerExpr={
    if(is.null(input$available_seurats_rows_selected))
      return(NULL)

    # row number is saved in the `input` so get the expression to `get` the object from the initial search table
    input_seurat_expr <- server_env$available_seurat_objects %>% pluck('choiceValue') %>% pluck(input$available_seurats_rows_selected)

    progress <- shiny::Progress$new(session=session, min=0, max=9/10)
    on.exit(progress$close())
    progress$set(value=0, message='Loading environment')

    # load Seurat object from user
    progress$inc(detail='Updating UI with cluster options')
    seurat <- parse(text=input_seurat_expr) %>% eval()
    seurat <- subset(seurat, subset=nFeature_RNA>0 & nCount_RNA>0)
    cluster_options <- c('seurat_clusters', str_subset(colnames(seurat@meta.data), '_snn_res.'))
    updateSelectInput(session=session, inputId='seurat_cluster_set.dd', choices=cluster_options)
    updateSelectInput(session=session, inputId='features_heatmap.seurat_cluster_set.dd', choices=cluster_options)
    for(nsid in module_environments$cluster_resolution_pickers$ns)
      updateSelectInput(session=session, inputId=nsid, choices=cluster_options)

    progress$inc(detail='Checking meta.data')
    if(is.null(seurat@meta.data$seurat_clusters))
      seurat@meta.data$seurat_clusters <- 0

    progress$inc(detail='Initialising reduced dimension plot')
    dimred_method <- Seurat:::DefaultDimReduc(seurat)
    if(!is.null(input$reduction_selection.dd))
      dimred_method <- input$reduction_selection.dd
    
    seurat@reductions[[dimred_method]]@cell.embeddings[,1:2] %>%
        as.data.frame() %>%
        set_names(c('DIMRED_1','DIMRED_2')) %>%
        cbind(seurat@meta.data) -> seurat_object.reactions$dimred

    progress$inc(detail='Counting clusters identified in each set')
    select_at(seurat@meta.data, vars(contains('_snn_res.'), 'seurat_clusters')) %>%
      mutate_all(function(x) {as.character(x) %>% as.numeric()}) %>%
      gather(key='cluster_set', value='ID') %>%
      group_by(cluster_set) %>%
      summarise(N=length(unique(ID))) %>%
      deframe() -> clusters_per_resolution

    progress$inc(detail='Getting summary statistics')
    list(n_cells=nrow(seurat@meta.data),
         total_reads=sum(seurat@meta.data$nCount_RNA),
         median_reads_per_cell=round(x=median(seurat@meta.data$nCount_RNA), digits=0),
         median_genes_per_cell=round(x=median(seurat@meta.data$nFeature_RNA), digits=0),
         min_reads_per_cell=min(seurat@meta.data$nCount_RNA), max_reads_per_cell=max(seurat@meta.data$nCount_RNA),
         min_genes_per_cell=min(seurat@meta.data$nFeature_RNA), max_genes_per_cell=max(seurat@meta.data$nFeature_RNA),
         max_percent_mitochondria=round(max(seurat@meta.data$percent_mt)+0.05, digits=1)) -> cell_filtering_data.reference

    progress$inc(detail='Setting default assay')
    selected_assay <- 'RNA'
    DefaultAssay(seurat) <- selected_assay
    if(sum(seurat@assays[[selected_assay]]@counts)==sum(seurat@assays[[selected_assay]]@data))
      seurat <- NormalizeData(seurat)

    progress$inc(detail='Updating UI elements')
    updateTextInput(session=session, inputId='min_features_per_cell.textinput', placeholder=cell_filtering_data.reference$min_genes_per_cell)
    updateTextInput(session=session, inputId='max_features_per_cell.textinput', placeholder=cell_filtering_data.reference$max_genes_per_cell)
    updateTextInput(session=session, inputId='min_expression_per_cell.textinput', placeholder=cell_filtering_data.reference$min_reads_per_cell)
    updateTextInput(session=session, inputId='max_expression_per_cell.textinput', placeholder=cell_filtering_data.reference$max_reads_per_cell)
    updateTextInput(session=session, inputId='percent_mitochondria.textinput', placeholder=cell_filtering_data.reference$max_percent_mitochondria)
    updateSelectInput(session=session, inputId='reduction_selection.dd', choices=names(seurat@reductions), selected=Seurat:::DefaultDimReduc(seurat))
    updateSelectInput(session=session, inputId='assay_selection.dd', choices=Assays(seurat), selected=DefaultAssay(seurat))
    update_autocomplete_input(session=session, id='gene_of_interest.dd', options=c(sort(rownames(seurat)), 'nFeature_RNA', 'nCount_RNA', 'percent_mt', 'orig.ident', 'orig.species','orig.timepoint','orig.tissue','orig.replicate'))

    progress$inc(detail='Saving variables')
    available_assays <- Assays(seurat)
    available_slots <- lapply(seurat@assays, function(x) c('counts','data','scale.data') %>% purrr::set_names() %>% lapply(function(y) slot(x,y) %>% nrow())) %>% lapply(function(y) names(y)[unlist(y)>0])

   # copy the important stuff into the reaction values
    seurat_object.reactions$active_object_expr <- input_seurat_expr
    seurat_object.reactions$seurat <- seurat
    seurat_object.reactions$mart <- seurat@misc$mart
    seurat_object.reactions$formatted.project.name <- seurat@project.name %>% str_replace_all(pattern='_', replacement=' ') %>% str_to_upper()
    seurat_object.reactions$reference_metrics <- cell_filtering_data.reference
    seurat_object.reactions$clusters_per_resolution <- clusters_per_resolution
    seurat_object.reactions$selected_clusters_per_resolution <- clusters_per_resolution[cluster_options[1]]
  })
}
