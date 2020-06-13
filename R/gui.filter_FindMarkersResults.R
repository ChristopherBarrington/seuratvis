
filter_FindMarkersResults.ui <- function(id) {
  sprintf(fmt='### %s-filter_FindMarkersResults.ui', id) %>% message()

  module <- 'filter_FindMarkersResults'

  # make unique id for this object
  ns <- NS(namespace=id)
  module_ns <- ns(id=module)

  # create an environment in the seuratvis namespace
  e <- new.env()
  e$id <- id
  assign(x=module_ns, val=e, envir=module_environments)

  # record the server(s) to call
  module_environments$filter_FindMarkersResults$id %<>% c(id) # keep track so other modules can update
  get0(env=module_servers_to_call, x=id) %>% append(sprintf(fmt='%s.server', module)) %>% assign(env=module_servers_to_call, x=id)

  # make ui elements
  filtering_sidebar <- filterDF_UI(id=ns('filter_parameters'))
  dimred_plot <- esquisserUI(id=ns('dimred_plot'), 
                             header=FALSE, choose_data=FALSE,
                             insert_code=FALSE, disable_filters=TRUE)

  # return ui element(s)
  list(fluidRow(column(width=2, filtering_sidebar),
                column(width=10, DT::dataTableOutput(outputId=ns('table')),
                       hr(),
                       column(width=9, dimred_plot))))
}

filter_FindMarkersResults.server <- function(input, output, session) {
  session$ns('') %>% sprintf(fmt='### %sfilter_FindMarkersResults.server') %>% message()

  # get environments containing variables to run/configure this object
  collect_environments(id=parent.frame()$id, module='filter_FindMarkersResults') # provides `seuratvis_env`, `server_env` and `module_env`
  session_server <- get(x='session', env=server_env)

  # update UI when Seurat object is loaded
  observeEvent(eventExpr=seurat_object.reactions$seurat, handlerExpr={
    # send a message
    session$ns('') %>% sprintf(fmt='### %sfilter_FindMarkersResults.server-observeEvent-seurat_object.reactions$seurat [%s]', seurat_object.reactions$formatted.project.name) %>% message()

    # create varaibles for shorthand
    seurat <- seurat_object.reactions$seurat
    values <- FetchData(seurat, c('UMAP_1', 'UMAP_2', 'tSNA_1', 'tSNA_2', str_subset(colnames(seurat@meta.data), '_snn_res.')))
    seurat@misc$FindMarkersResults %>%
      pluck('wilcox') %>%
      mutate(p_adj_group={p_val_adj %>% cut(breaks=c(0, 0.1/100, 1/100, 5/100, 10/100, 100/100), labels=c('<0.1%','<1%','<5%','<10%','NS'), include.lowest=TRUE, right=TRUE)}) %>%
      dplyr::select(cluster_set, ident.1, gene, pct.1, pct.2, avg_logFC, p_adj_group) %>%
      rename(`Cluster set`='cluster_set', `Cluster ID`='ident.1', `Gene`='gene', `Cluster detection`='pct.1', `Map detection`='pct.2', `Avg. logFC`='avg_logFC', `Adj. P`='p_adj_group') -> tidied_results

    # update the reactive(s)
    seurat_object.reactions$FindMarkersResults$table <- tidied_results
    seurat_object.reactions$FindMarkersResults$vars <- c('Cluster set', 'Cluster ID', 'Adj. P', 'Avg. logFC', 'Cluster detection', 'Map detection')

    # call the ggplot module for the reduced dimension data
    callModule(session=session_server,
               module=esquisserServer,
               id=session$ns('dimred_plot'),
               data=reactiveValues(data=values, name=seurat_object.reactions$formatted.project.name))})

  # handle the data.frame filtering
  ## call the data.frame filtering filtering module
  callModule(session=session_server, module=filterDF, id=session$ns('filter_parameters'),
             picker=TRUE,
             data_name=reactive(seurat_object.reactions$formatted.project.name), 
             data_table=reactive(seurat_object.reactions$FindMarkersResults$table),
             data_vars=reactive(seurat_object.reactions$FindMarkersResults$vars)) -> res_filter
  
  ## render the filtered data.table
  DT::renderDataTable({
    # send a message
    session$ns('') %>% sprintf(fmt='### %sfilter_FindMarkersResults.server-DT::renderDataTable [%s]', seurat_object.reactions$formatted.project.name) %>% message()

    # render the datatable and send it to the output
    res_filter$data_filtered() %>%
      DT::datatable(rownames=FALSE,
                    options=list(columnDefs=list(list(className='dt-right', targets=c(1,6))),
                                 ordering=FALSE,
                                 dom='litp'),
                    style='bootstrap4',
                    class='stripe') %>%
      formatRound(columns=c('Cluster detection', 'Map detection'), digits=2) %>%
      formatRound(columns=c('Avg. logFC'), digits=3)}) -> output$table
}
