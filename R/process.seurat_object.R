#' Load a seurat object into a reactive
#' 
#' @rdname available_seurats
#' 
seurat_object.server <- function(input, output, session, seurat, ...) {
  session$ns('') %>% sprintf(fmt='### %sload_a_seurat.server') %>% message()

  # get environments containing variables to run/configure this object
  collect_environments(id=parent.frame()$id, module='available_seurats') # provides `seuratvis_env`, `server_env` and `module_env`
  session_server <- get(x='session', env=server_env)

  # make a seurat object reactives list to keep all of the information
  seurat.rv <- reactiveValues()

  # react to a row being selected in the available Seurats table
  observeEvent(eventExpr=input$seurats_table_rows_selected, handlerExpr={
    # make sure these elements are defined
    req(input$seurats_table_rows_selected)

    # send a message
    session$ns('') %>% sprintf(fmt='### %sload_a_seurat.server-observeEvent-input$seurats_table_rows_selected [%s]', input$seurats_table_rows_selected) %>% message()

    # clear out seurat.rv
    #! TODO: better way to reset the reactive values list?
    for(i in names(seurat.rv)) seurat.rv[[i]] <- NULL

    # row number is saved in the `input` so get the expression to `get` the object from the initial search table
    input_seurat_expr <- server_env$available_seurat_objects %>% pluck('choiceValue') %>% pluck(input$seurats_table_rows_selected)

    # load Seurat object from user
    seurat <- parse(text=input_seurat_expr) %>% eval()

    # set some dummy clusters if there are none
    if(is.null(seurat@meta.data$seurat_clusters))
      seurat@meta.data$seurat_clusters <- 0

    # ensure we are using RNA assay and it is normalised
    selected_assay <- 'RNA'
    DefaultAssay(seurat) <- selected_assay
    if(sum(seurat@assays[[selected_assay]]@counts)==sum(seurat@assays[[selected_assay]]@data))
      seurat <- NormalizeData(seurat)

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

    # pull gene modules out of the meta.data
    meta_data <- gene_modules <- seurat@meta.data
    gene_modules_regex <- '^GeneModule-'
    gene_modules %<>% (function(x) x[colnames(x) %>% str_detect(regex(pattern=gene_modules_regex, ignore_case=TRUE), negate=FALSE)]) %>% set_names(str_remove, gene_modules_regex)
    meta_data %<>% (function(x) x[colnames(x) %>% str_detect(regex(pattern=gene_modules_regex, ignore_case=TRUE), negate=TRUE)])

    # copy the important stuff into the reaction values
    seurat.rv$active_object_expr <- input_seurat_expr
    seurat.rv$mart <- Misc(object=seurat, slot='mart')
    seurat.rv$n_cells <- ncol(seurat)
    seurat.rv$n_features <- nrow(seurat) #! TODO: this is the number of features in the _active_ assay
    seurat.rv$metadata <- meta_data
    seurat.rv$gene_modules <- gene_modules
    seurat.rv$project <- Project(seurat)
    seurat.rv$formatted_project <- Project(seurat) %>% str_replace_all(pattern='_', replacement=' ') %>% str_to_upper()
    seurat.rv$object <- seurat

    seurat.rv$done <- rnorm(1)})

  # react when the configuration options are changed
  ## react to the percent mitochondria column being set
  observeEvent(eventExpr=c(input$proportion_mt_picker, seurat.rv$object), handlerExpr={
    # make sure these elements are defined
    req(seurat.rv$object)
    req(input$proportion_mt_picker)

    # send a message
    session$ns('') %>% sprintf(fmt='### %sload_a_seurat.server-observeEvent-input$proportion_mt_picker [%s]', input$proportion_mt_picker) %>% message()

    # create varaibles for shorthand
    seurat <- seurat.rv$object
    var <- input$proportion_mt_picker
    values <- FetchData(object=seurat, vars=var) %>% set_names('proportion_mt')
    high <- max(values$proportion_mt)

    # trigger update of the ui element(s)
    seurat.rv$reset_proportion_mt <- rnorm(1)

    # update the reactive
    seurat.rv$proportion_mt_values_max <- high
    seurat.rv$proportion_mt_values <- values
    seurat.rv$proportion_mt_variable <- var})
  
  ## react to the number of features per cell column being set
  observeEvent(eventExpr=c(input$n_features_picker, seurat.rv$object), handlerExpr={
    # make sure these elements are defined
    req(seurat.rv$object)
    req(input$n_features_picker)

    # send a message
    session$ns('') %>% sprintf(fmt='### %sload_a_seurat.server-observeEvent-input$n_features_picker [%s]', input$n_features_picker) %>% message()

    # create varaibles for shorthand
    seurat <- seurat.rv$object
    var <- input$n_features_picker
    values <- FetchData(object=seurat, vars=var) %>% set_names('n_features')
    low <- min(values)
    high <- max(values)

    # trigger update of the ui element(s)
    seurat.rv$reset_n_features <- rnorm(1)

    # update the reactive
    seurat.rv$n_features_values_min <- low
    seurat.rv$n_features_values_max <- high
    seurat.rv$n_features_values <- values
    seurat.rv$n_features_variable <- var})

  ## react to the number of UMI column being set
  observeEvent(eventExpr=c(input$n_umi_picker, seurat.rv$object), handlerExpr={
    # make sure these elements are defined
    req(seurat.rv$object)
    req(input$n_umi_picker)

    # send a message
    session$ns('') %>% sprintf(fmt='### %sload_a_seurat.server-observeEvent-input$n_umi_picker [%s]', input$n_umi_picker) %>% message()

    # create varaibles for shorthand
    seurat <- seurat.rv$object
    var <- input$n_umi_picker
    values <- FetchData(object=seurat, vars=var) %>% set_names('n_umi')
    low <- min(values)
    high <- max(values)

    # trigger update of the ui element(s)
    seurat.rv$reset_n_umi <- rnorm(1)

    # update the reactives
    seurat.rv$n_umi <- sum(values)
    seurat.rv$n_umi_values_min <- low
    seurat.rv$n_umi_values_max <- high
    seurat.rv$n_umi_values <- values
    seurat.rv$n_umi_variable <- var})
 
  # return the reactive values list
  return(seurat.rv)
}
