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
  # get0(env=module_servers_to_call, x=id) %>% append(sprintf(fmt='%s.server', c('available_seurats', 'load_a_seurat'))) %>% assign(env=module_servers_to_call, x=id)

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
available_seurats.server <- function(input, output, session, ...) {
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
    `Reductions`='reductions',
    `Size`='size') -> column_order

  # make the `data.frame` of Seurat information
  collapse_strings <- function(x, replacement='-', sep=', ')
    ifelse(length(x)==0 | (length(x)==1 && x==''), replacement, str_c(x, collapse=sep))

  if(nrow(server_env$available_seurat_objects)==0) {
    #! TODO: use confirmSweetAlert to send a message to an observer to `stop('no seurats loaded!')` the app
    sendSweetAlert(session=session, type='error', html=FALSE, title='We have a problem!', text='No Seurat objects found!', btn_labels='OK')
    return(NULL)
  }

  # make the data.frame for the detected Seurat objects
  server_env$available_seurat_objects %>%
    dplyr::select(-choiceName) %>%
    plyr::adply(.margins=1, .parallel=FALSE, function(params) { #! TODO: parallelise this doMC does not work on Windows though!
      x <- eval(parse(text=params$choiceValue))
      
      # choose a gene to look for to identify male dataset
      x@misc$mart@dataset %>%
        str_remove(pattern='_.*') %>%
        switch(hsapiens='SRY',
               mmusculus='Sry',
               'SRY') -> boy_gene

      # make the data.frame for this Seurat
      data.frame(project=reformat_project_name(Project(x)),
                 environment=if_else(params$env=='globalenv()', 'RGlobal', as.character(params$env)),
                 ncells=ncol(x),
                 numi=sum(x@meta.data$nCount_RNA),
                 median_umi=median(x@meta.data$nCount_RNA),
                 dimensions=ifelse(is.null(x@misc$n_dimensions$pca), -1, x@misc$n_dimensions$pca),
                 filtered=!is.null(x@misc$cells_filtered) && x@misc$cells_filtered,
                 integrated=!is.null(x@misc$integrated_dataset) && x@misc$integrated_dataset,
                 active_assay=DefaultAssay(x),
                 assays={Assays(x) %>% str_subset(pattern=DefaultAssay(x), negate=TRUE) %>% collapse_strings()},
                 nfeatures=nrow(x),
                 reductions={Reductions(x) %>% collapse_strings()},
                 guessed_sex={FetchData(x, vars=boy_gene) %>% is_greater_than(0) %>% any() %>% if_else(as.character(icon(name='mars', class='boy')), as.character(icon(name='venus', class='girl')))},
                 size={object.size(x) %>% format(units='Gb') %>% str_remove(' ')})}) %>%
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
      formatStyle(columns=c('object', 'size'),
                  fontFamily='monospace',
                  fontWeight='bold') %>%
      formatRound(columns=c('ncells', 'numi', 'median_umi', 'nfeatures'),
                  digits=0) %>%
      formatStyle(columns={sapply(data_to_show, class) %>% str_which('logical')},
                  color=styleEqual(levels=c(0,1), values=c('#B96A3C','#B93C8B'), default='orange'),
                  fontFamily='monospace',
                  fontWeight='bold')}) -> output$seurats_table
}
