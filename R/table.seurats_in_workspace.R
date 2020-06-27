#'
#' 
seurats_in_workspace.table <- function(id) {
  DTOutput(outputId=NS(id, 'table')) %>% withSpinner(type=6)
}

#'
#' 
seurats_in_workspace.server <- function(input, output, session) {
  # look in the workspace for Seurat objects
  find_seurat_objects <- function() {
    ls(envir=globalenv()) %>%
      sapply(function(O) get(x=O, envir=globalenv()) %>% class()) %>%
      unlist() %>%
      enframe() %>%
      plyr::dlply(~value, pluck, 'name') %>%
      pluck('environment') %>% # select the environments from the RGlobalEnv
      rev() %>%
      sapply(get, envir=globalenv()) %>%
      append(list(`globalenv()`=globalenv())) %>% # make a list of all environments and RGlobalEnv
      rev() %>%
      lapply(function(E) {ls(envir=E) %>% sapply(function(O) get(x=O, envir=E) %>% class()) %>% unlist() %>% enframe() %>% plyr::dlply(~value, pluck, 'name') %>% pluck('Seurat')}) %>% # get the class of objects in the environments and select out the Seurat objects
      plyr::ldply(.id='env', enframe) %>%
      dplyr::select(-name) %>%
      unite(col='choiceValue', sep='$', env, value, remove=FALSE) %>%
      (function(x) {
        if(length(unique(x$env))==1) {
          x %>% mutate(choiceName=value)
        } else {
          x %>% mutate(choiceName=sprintf('%s [%s]', str_replace_all(string=value, pattern='_', replacement=' '), str_remove_all(string=str_replace_all(string=env, pattern='_', replacement=' '), pattern='\\(\\)$')))
        }}) %>%
      arrange(choiceName) %>%
      mutate(env=as.character(env)) -> available_objects # add varaibles for where to find the objects and what to call them
  }

  ## run the search function
  find_seurat_objects() -> available_seurat_objects

  ## check that there are some Seurats in the workspace
  if(nrow(available_seurat_objects)==0) {
    #! TODO: use confirmSweetAlert to send a message to an observer to `stop('no seurats loaded!')` the app
    sendSweetAlert(session=session, type='error', html=FALSE, title='We have a problem!', text='No Seurat objects found!', btn_labels='OK')
    return(NULL)
  }

  # make the `data.frame` of Seurat information
  ## this function will join strings or retun '-' if there is nothing to join
  collapse_strings <- function(x, replacement='-', sep=', ')
    ifelse(length(x)==0 | (length(x)==1 && x==''), replacement, str_c(x, collapse=sep))

  ## define output column names and order
  c(Project='project',
    Environment='environment',
    `Object name`='object',
    # `Size`='size',
    # `Guessed sex`='guessed_sex',
    `Number of cells`='ncells',
    `Total UMI`='numi',
    `Median UMI`='median_umi',
    `Filtered?`='filtered',
    `Integrated?`='integrated',
    `Active assay`='active_assay',
    `Other assays`='assays',
    # `Features in active assay`='nfeatures',
    `Reductions`='reductions') -> column_order

  ## make the data.frame for the detected Seurat objects
  available_seurat_objects %>%
    dplyr::select(-choiceName) %>%
    plyr::adply(.margins=1, .parallel=FALSE, function(params) { #! TODO: parallelise this doMC does not work on Windows though!
      # make a copy of the Seurat object
      x <- eval(parse(text=params$choiceValue))

      # make the data.frame for this Seurat
      data.frame(project=reformat_project_name(Project(x)),
                 environment=if_else(params$env=='globalenv()', 'RGlobal', as.character(params$env)),
                 ncells=ncol(x),
                 numi=sum(x@meta.data$nCount_RNA),
                 median_umi=median(x@meta.data$nCount_RNA),
                 filtered=!is.null(x@misc$cells_filtered) && x@misc$cells_filtered,
                 integrated=!is.null(x@misc$integrated_dataset) && x@misc$integrated_dataset,
                 active_assay=DefaultAssay(x),
                 assays={Assays(x) %>% str_subset(pattern=DefaultAssay(x), negate=TRUE) %>% collapse_strings()},
                 # nfeatures=nrow(x),
                 reductions={Reductions(x) %>% collapse_strings()})}) %>%
    mutate(object=value) %>%
    select_at(vars(all_of(column_order), everything())) -> data_to_show

  # make and format the `datatable`
  ## identify columns to format in the `datatable` `columnDefs` argument
  formatted_colnames <- c(names(column_order), colnames(data_to_show)[! colnames(data_to_show) %in% column_order])
  list(hide={seq(length(column_order)+1, ncol(data_to_show))},
       center={c(sapply(data_to_show, class) %>% str_which('logical'),
                 str_which(column_order, 'guessed_sex'))}) %>%
    lapply(subtract, e2=1) -> columnDef_targets

  ## render the datatable
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
                    selection=list(mode='none'),
                    escape=FALSE) %>%
      formatStyle(columns='ncells',
                  background=styleColorBar(data=c(0,max(data_to_show$ncells)), color='#605CA8'),
                  backgroundSize='98% 50%',
                  backgroundRepeat='no-repeat',
                  backgroundPosition='center') %>%
      formatStyle(columns='numi',
                  background=styleColorBar(data=c(0,max(data_to_show$numi)), color='#605CA8'),
                  backgroundSize='98% 50%',
                  backgroundRepeat='no-repeat',
                  backgroundPosition='center') %>%
      formatStyle(columns='median_umi',
                  background=styleColorBar(data=c(0,max(data_to_show$median_umi)), color='#605CA8'),
                  backgroundSize='98% 50%',
                  backgroundRepeat='no-repeat',
                  backgroundPosition='center') %>%
      # formatStyle(columns='nfeatures',
      #             background=styleColorBar(data=c(0,max(data_to_show$nfeatures)), color='#3CB96A'),
      #             backgroundSize='98% 50%',
      #             backgroundRepeat='no-repeat',
      #             backgroundPosition='center') %>%
      formatStyle(columns=c('object'),
                  fontFamily='monospace',
                  fontWeight='bold') %>%
      formatRound(columns=c('ncells', 'numi', 'median_umi'),
                  digits=0) %>%
      formatStyle(columns={sapply(data_to_show, class) %>% str_which('logical')},
                  color=styleEqual(levels=c(0,1), values=c('#dd4b39','#00A65A'), default='orange'),
                  fontFamily='monospace',
                  fontWeight='bold')}) -> output$table

  # return the table of available Seurat objects
  available_seurat_objects
}