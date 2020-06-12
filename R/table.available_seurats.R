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
#' @import DT
#' 
#' @rdname available_seurats
#' 
available_seurats.ui <- function(id) {
  sprintf(fmt='### %s-available_seurats.ui', id) %>% message()

  module <- 'available_seurats'

  # make unique id for this object
  ns <- NS(namespace=id)
  module_ns <- ns(id=module)

  # create an environment in the seuratvis namespace
  e <- new.env()
  e$id <- id
  assign(x=module_ns, val=e, envir=module_environments)

  # record the server(s) to call
  get0(env=module_servers_to_call, x=id) %>% append(sprintf(fmt='%s.server', c('available_seurats', 'load_a_seurat'))) %>% assign(env=module_servers_to_call, x=id)

  # make ui elements to pick columns to use
  selectizeInput_defaults <- list(choices=NULL, selected=NULL, multiple=FALSE)
  list(inputId=ns(id='n_features_picker'), label='Features') %>% modifyList(x=selectizeInput_defaults) %>% do.call(what=selectizeInput) -> n_features_picker
  list(inputId=ns(id='n_umi_picker'), label='UMI') %>% modifyList(x=selectizeInput_defaults) %>% do.call(what=selectizeInput) -> n_umi_picker
  list(inputId=ns(id='proportion_mt_picker'), label='Mitochondrial proportion') %>% modifyList(x=selectizeInput_defaults) %>% do.call(what=selectizeInput) -> proportion_mt_picker

  # return ui element(s)
  tagList(DTOutput(outputId=ns(id='seurats_table')), tags$h3('Configure metadata columns'),
          column(width=6, splitLayout(n_features_picker, n_umi_picker, proportion_mt_picker),
          tags$head(tags$style(HTML(".shiny-split-layout > div {overflow: visible; }")))))
}

#' Render the table of Seurat object and their metadata
#' 
#' @import DT
#' 
#' @rdname available_seurats
#' 
available_seurats.server <- function(input, output, session) {
  session$ns('') %>% sprintf(fmt='### %savailable_seurats.server') %>% message()

  # get environments containing variables to run/configure this object
  collect_environments(id=parent.frame()$id, module='available_seurats') # provides `seuratvis_env`, `server_env` and `module_env`

  # define output column names and order
  c(Project='project',
    Environment='environment',
    `Object name`='object',
    `Guessed sex`='guessed_sex',
    `Number of cells`='ncells',
    `Total UMI`='numi',
    `Median UMI`='median_umi',
    `Filtered?`='filtered',
    `Integrated?`='integrated',
    `Active assay`='active_assay',
    `Other assays`='assays',
    `Features in active assay`='nfeatures',
    `Reductions`='reductions') -> column_order

  # make the `data.frame` of Seurat information
  collapse_strings <- function(x, replacement='-', sep=', ')
    ifelse(length(x)==0 | (length(x)==1 && x==''), replacement, str_c(x, collapse=sep))

  if(nrow(server_env$available_seurat_objects)==0) {
    #! TODO: use confirmSweetAlert to send a message to an observer to `stop('no seurats loaded!')` the app
    sendSweetAlert(session=session, type='error', html=FALSE, title='We have a problem!', text='No Seurat objects found!', btn_labels='OK')
    return(NULL)
  }

  server_env$available_seurat_objects %>%
    dplyr::select(-choiceName) %>%
    plyr::adply(.margins=1, function(params) {
      x <- eval(parse(text=params$choiceValue))
      data.frame(project=reformat_project_name(Project(x)),
                 environment=if_else(params$env=='globalenv()', 'RGlobal', as.character(params$env)),
                 ncells=ncol(x),
                 numi=sum(x@meta.data$nCount_RNA),
                 median_umi=median(x@meta.data$nCount_RNA),
                 dimensions=ifelse(is.null(x@misc$n_dimensions), -1, x@misc$n_dimensions),
                 filtered=!is.null(x@misc$cells_filtered) && x@misc$cells_filtered,
                 integrated=!is.null(x@misc$integrated_dataset) && x@misc$integrated_dataset,
                 active_assay=DefaultAssay(x),
                 assays={Assays(x) %>% str_subset(pattern=DefaultAssay(x), negate=TRUE) %>% collapse_strings()},
                 nfeatures=nrow(x),
                 reductions={Reductions(x) %>% collapse_strings()},
                 guessed_sex={FetchData(x, vars=c('SRY','Sry')) %>% is_greater_than(0) %>% any() %>% if_else(as.character(icon(name='mars', class='boy')), as.character(icon(name='venus', class='girl')))})}) %>%
    mutate(dimensions=as.integer(dimensions),
           object=value) %>%
    select_at(vars(all_of(column_order), everything())) -> data_to_show

  # identify columns to format in the `datatable` `columnDefs` argument
  formatted_colnames <- c(names(column_order), colnames(data_to_show)[! colnames(data_to_show) %in% column_order])
  list(hide={seq(length(column_order)+1, ncol(data_to_show))},
       center={c(sapply(data_to_show, class) %>% str_which('logical'),
                 str_which(column_order, 'guessed_sex'))}) %>%
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
                                   selected=input$seurats_table_rows_selected),
                    escape=FALSE) %>%
      formatStyle(columns='ncells',
                  background=styleColorBar(data=c(0,max(data_to_show$ncells)), color='#3CB96A'),
                  backgroundSize='98% 50%',
                  backgroundRepeat='no-repeat',
                  backgroundPosition='center') %>%
      formatStyle(columns='numi',
                  background=styleColorBar(data=c(0,max(data_to_show$numi)), color='#3CB96A'),
                  backgroundSize='98% 50%',
                  backgroundRepeat='no-repeat',
                  backgroundPosition='center') %>%
      formatStyle(columns='median_umi',
                  background=styleColorBar(data=c(0,max(data_to_show$median_umi)), color='#3CB96A'),
                  backgroundSize='98% 50%',
                  backgroundRepeat='no-repeat',
                  backgroundPosition='center') %>%
      formatStyle(columns='nfeatures',
                  background=styleColorBar(data=c(0,max(data_to_show$nfeatures)), color='#3CB96A'),
                  backgroundSize='98% 50%',
                  backgroundRepeat='no-repeat',
                  backgroundPosition='center') %>%
      formatStyle(columns='object',
                  fontFamily='monospace',
                  fontWeight='bold') %>%
      formatRound(columns=c('ncells', 'numi', 'median_umi', 'nfeatures'),
                  digits=0) %>%
      formatStyle(columns={sapply(data_to_show, class) %>% str_which('logical')},
                  color=styleEqual(levels=c(0,1), values=c('#B96A3C','#B93C8B'), default='orange'),
                  fontFamily='monospace',
                  fontWeight='bold')}) -> output$seurats_table
}

#' Load a seurat object into a reactive
#' 
#' @rdname available_seurats
#' 
load_a_seurat.server <- function(input, output, session) {
  session$ns('') %>% sprintf(fmt='### %sload_a_seurat.server') %>% message()

  # get environments containing variables to run/configure this object
  collect_environments(id=parent.frame()$id, module='available_seurats') # provides `seuratvis_env`, `server_env` and `module_env`
  session_server <- get(x='session', env=server_env)

  # react to a row being selected in the available Seurats table
  observeEvent(eventExpr=input$seurats_table_rows_selected, handlerExpr={
    # make sure these elements are defined
    req(input$seurats_table_rows_selected)

    # send a message
    session$ns('') %>% sprintf(fmt='### %sload_a_seurat.server-observeEvent-input$seurats_table_rows_selected [%s]', input$seurats_table_rows_selected) %>% message('')

    # empty out the reactives
    for(i in names(seurat_object.reactions))
      seurat_object.reactions[[i]] <- NULL

    for(i in names(filtering_parameters.reactions))
      filtering_parameters.reactions[[i]] <- NULL

    for(i in names(filtered_cells.reactions))
      filtered_cells.reactions[[i]] <- NULL

    for(i in names(seurat_configuration.reactions))
      seurat_configuration.reactions[[i]] <- NULL

    for(i in names(filtering_arguments.reactions))
      filtering_arguments.reactions[[i]] <- NULL

    for(i in names(reference_metrics.rv))
      reference_metrics.rv[[i]] <- NULL

    for(i in names(plotting_options.rv))
      plotting_options.rv[[i]] <- NULL

    for(i in names(selections.rv))
      selections.rv[[i]] <- NULL

    # row number is saved in the `input` so get the expression to `get` the object from the initial search table
    input_seurat_expr <- server_env$available_seurat_objects %>% pluck('choiceValue') %>% pluck(input$seurats_table_rows_selected)

    # load Seurat object from user
    seurat <- parse(text=input_seurat_expr) %>% eval()

    if(is.null(seurat@meta.data$seurat_clusters))
      seurat@meta.data$seurat_clusters <- 0

    selected_assay <- 'RNA'
    DefaultAssay(seurat) <- selected_assay
    if(sum(seurat@assays[[selected_assay]]@counts)==sum(seurat@assays[[selected_assay]]@data))
      seurat <- NormalizeData(seurat)

    #! TODO: remove seurat_object.reactions$reference_metrics and cell_filtering_data.reference
    list(n_cells=ncol(seurat),
         n_features=nrow(seurat),
         # total_umi=sum(seurat@meta.data$nCount_RNA),
         # median_umi_per_cell=round(x=median(seurat@meta.data$nCount_RNA), digits=0),
         # median_features_per_cell=round(x=median(seurat@meta.data$nFeature_RNA), digits=0),
         # min_umi_per_cell=min(seurat@meta.data$nCount_RNA), max_umi_per_cell=max(seurat@meta.data$nCount_RNA),
         # min_features_per_cell=min(seurat@meta.data$nFeature_RNA), max_features_per_cell=max(seurat@meta.data$nFeature_RNA),
         
         # total_umi_per_cell_min=min(seurat@meta.data$nCount_RNA), total_umi_per_cell_max=max(seurat@meta.data$nCount_RNA),
         # features_per_cell_min=min(seurat@meta.data$nFeature_RNA), features_per_cell_max=max(seurat@meta.data$nFeature_RNA),
         project=Project(seurat)) -> cell_filtering_data.reference

    available_assays <- Assays(seurat)
    available_slots <- lapply(seurat@assays, function(x) c('counts','data','scale.data') %>% purrr::set_names() %>% lapply(function(y) slot(x,y) %>% nrow())) %>% lapply(function(y) names(y)[unlist(y)>0])

    # update ui elements
    ## get the numeric metadata variables
    sapply(seurat@meta.data, is.numeric) %>% subset(x=colnames(seurat@meta.data)) -> choices

    ## guess a default choice
    n_features_picker_default <- preferred_choice(x=choices, preferences=c('nFeature_RNA','nFeature_SCT'))
    n_umi_picker_default <- preferred_choice(x=choices, preferences=c('nCount_RNA','nCount_SCT'))
    proportion_mt_picker_default <- preferred_choice(x=choices, preferences=c('percent.mt', 'percent_mt', 'prop.mt', 'prop_mt'))

    ## define the choices and default in the input ui elements
    updateSelectizeInput(session=session, inputId='n_features_picker', choices=choices, selected=n_features_picker_default)
    updateSelectizeInput(session=session, inputId='n_umi_picker', choices=choices, selected=n_umi_picker_default)
    updateSelectizeInput(session=session, inputId='proportion_mt_picker', choices=choices, selected=proportion_mt_picker_default)

    # copy the important stuff into the reaction values
    seurat_object.reactions$active_object_expr <- input_seurat_expr
    seurat_object.reactions$mart <- seurat@misc$mart
    seurat_object.reactions$formatted.project.name <- seurat@project.name %>% str_replace_all(pattern='_', replacement=' ') %>% str_to_upper()
    seurat_object.reactions$reference_metrics <- cell_filtering_data.reference
    seurat_object.reactions$cell_metadata <- seurat@meta.data
    seurat_object.reactions$project <- Project(seurat)
    seurat_object.reactions$seurat <- seurat})

  # react when the configuration options are changed
  ## react to the percent mitochondria column being set
  observeEvent(eventExpr=c(input$proportion_mt_picker, seurat_object.reactions$seurat), handlerExpr={
    # make sure these elements are defined
    req(seurat_object.reactions$seurat)
    req(input$proportion_mt_picker)

    # send a message
    session$ns('') %>% sprintf(fmt='### %sload_a_seurat.server-observeEvent-input$proportion_mt_picker [%s]', input$proportion_mt_picker) %>% message('')

    # create varaibles for shorthand
    seurat <- seurat_object.reactions$seurat
    var <- input$proportion_mt_picker
    values <- FetchData(object=seurat, vars=var) %>% set_names('proportion_mt')
    high <- max(values$proportion_mt)

    # trigger update of the ui element(s)
    seurat_configuration.reactions$reset_proportion_mt <- rnorm(1)

    # update the reactive
    seurat_object.reactions$proportion_mt_values_max <- high
    seurat_object.reactions$proportion_mt_values <- values
    seurat_configuration.reactions$proportion_mt_variable <- var
    filtering_parameters.reactions$max_percent_mitochondria <- high})
  
  ## react to the number of features per cell column being set
  observeEvent(eventExpr=c(input$n_features_picker, seurat_object.reactions$seurat), handlerExpr={
    # make sure these elements are defined
    req(seurat_object.reactions$seurat)
    req(input$n_features_picker)

    # send a message
    session$ns('') %>% sprintf(fmt='### %sload_a_seurat.server-observeEvent-input$n_features_picker [%s]', input$n_features_picker) %>% message('')

    # create varaibles for shorthand
    seurat <- seurat_object.reactions$seurat
    var <- input$n_features_picker
    values <- FetchData(object=seurat, vars=var) %>% set_names('n_features')
    low <- min(values)
    high <- max(values)

    # trigger update of the ui element(s)
    seurat_configuration.reactions$reset_n_features <- rnorm(1)

    # update the reactive
    seurat_object.reactions$n_features_values_min <- low
    seurat_object.reactions$n_features_values_max <- high
    seurat_object.reactions$n_features_values <- values
    seurat_configuration.reactions$n_features_variable <- var
    filtering_parameters.reactions$features_per_cell_min <- low
    filtering_parameters.reactions$features_per_cell_max <- high})

  ## react to the number of UMI column being set
  observeEvent(eventExpr=c(input$n_umi_picker, seurat_object.reactions$seurat), handlerExpr={
    # make sure these elements are defined
    req(seurat_object.reactions$seurat)
    req(input$n_umi_picker)

    # send a message
    session$ns('') %>% sprintf(fmt='### %sload_a_seurat.server-observeEvent-input$n_umi_picker [%s]', input$n_umi_picker) %>% message('')

    # create varaibles for shorthand
    seurat <- seurat_object.reactions$seurat
    var <- input$n_umi_picker
    values <- FetchData(object=seurat, vars=var) %>% set_names('n_umi')
    low <- min(values)
    high <- max(values)

    # trigger update of the ui element(s)
    seurat_configuration.reactions$reset_n_umi <- rnorm(1)

    # update the reactives
    seurat_object.reactions$n_umi <- sum(values)
    seurat_object.reactions$n_umi_values_min <- low
    seurat_object.reactions$n_umi_values_max <- high
    seurat_object.reactions$n_umi_values <- values
    seurat_configuration.reactions$n_umi_variable <- var
    filtering_parameters.reactions$total_umi_per_cell_min <- low
    filtering_parameters.reactions$total_umi_per_cell_max <- high})
}
