#'
#' 
genes_in_modules.table <- function(id)
  DTOutput(outputId=NS(id,'genes_in_module'))

#' 
#' @import htmlwidgets
#' 
genes_in_modules.server <- function(input, output, session, seurat, picked_feature) {

  DT::renderDataTable({
  
    # make the data.frame to display
    ## get module(s)
    gm <- str_remove(picked_feature$name, 'GeneModule-') %>% str_split(pattern=',') %>% unlist()

    ## gene genes
    seurat$gene_modules[gm] %>%
      lapply(function(x) data.frame(`Gene.name`=as.character(x), stringsAsFactors=FALSE)) %>%
      lapply(purrr::transpose) -> genes_in_modules_list

    ## get a data.frame of gene module names
    seurat$gene_modules[gm] %>%
      names() %>%
      data.frame() %>%
      set_names('Gene module name') -> gene_module_names

    ## combine tha data.frame and list
    data <- cbind(` `='&oplus;', gene_module_names, `_details`=I(genes_in_modules_list))

    # define the callback
    #! TODO: look at this and figure it out!
    #! TODO: I don't like this implementation
    htmlwidgets::JS("table.column(1).nodes().to$().css({cursor: 'pointer'});",
       "",
       "// make the table header of the nested table",
       "var format = function(d, childId){",
       "  if(d != null){",
       "    var html = ", 
       "      '<table class=\"display compact hover\" id=\"' + childId + '\"><thead><tr>';",
       "    for (var key in d[d.length-1][0]) {",
       "      html += '<th>' + key + '</th>';",
       "    }",
       "    html += '</tr></thead></table>'",
       "    return html;",
       "  } else {",
       "    return '';",
       "  }",
       "};",
       "",
       "// row callback to style the rows of the child tables",
       "var rowCallback = function(row, dat, displayNum, index){",
       "  if($(row).hasClass('odd')){",
       "    $(row).css('background-color', 'papayawhip');",
       "    $(row).hover(function(){",
       "      $(this).css('background-color', '#E6FF99');",
       "    }, function() {",
       "      $(this).css('background-color', 'papayawhip');",
       "    });",
       "  } else {",
       "    $(row).css('background-color', 'lemonchiffon');",
       "    $(row).hover(function(){",
       "      $(this).css('background-color', '#DDFF75');",
       "    }, function() {",
       "      $(this).css('background-color', 'lemonchiffon');",
       "    });",
       "  }",
       "};",
       "",
       "// header callback to style the header of the child tables",
       "var headerCallback = function(thead, data, start, end, display){",
       "  $('th', thead).css({",
       "    'border-top': '3px solid indigo',", 
       "    'color': 'indigo',",
       "    'background-color': '#fadadd'",
       "  });",
       "};",
       "",
       "// make the datatable",
       "var format_datatable = function(d, childId){",
       "  var dataset = [];",
       "  var n = d.length - 1;",
       "  for(var i = 0; i < d[n].length; i++){",
       "    var datarow = $.map(d[n][i], function (value, index) {",
       "      return [value];",
       "    });",
       "    dataset.push(datarow);",
       "  }",
       "  var id = 'table#' + childId;",
       "  if (Object.keys(d[n][0]).indexOf('_details') === -1) {",
       "    var subtable = $(id).DataTable({",
       "                 'data': dataset,",
       "                 'autoWidth': true,",
       "                 'deferRender': true,",
       "                 'info': false,",
       "                 'lengthChange': false,",
       "                 'ordering': d[n].length > 1,",
       "                 'order': [],",
       "                 'paging': false,",
       "                 'scrollX': false,",
       "                 'scrollY': false,",
       "                 'searching': false,",
       "                 'sortClasses': false,",
       # "                 'rowCallback': rowCallback,",
       # "                 'headerCallback': headerCallback,",
       "                 'columnDefs': [{targets: '_all', className: 'dt-center'}]",
       "               });",
       "  } else {",
       "    var subtable = $(id).DataTable({",
       "            'data': dataset,",
       "            'autoWidth': true,",
       "            'deferRender': true,",
       "            'info': false,",
       "            'lengthChange': false,",
       "            'ordering': d[n].length > 1,",
       "            'order': [],",
       "            'paging': false,",
       "            'scrollX': false,",
       "            'scrollY': false,",
       "            'searching': false,",
       "            'sortClasses': false,",
       # "            'rowCallback': rowCallback,",
       # "            'headerCallback': headerCallback,",
       "            'columnDefs': [", 
       "              {targets: -1, visible: false},", 
       "              {targets: 0, orderable: false, className: 'details-control'},", 
       "              {targets: '_all', className: 'dt-center'}",
       "             ]",
       "          }).column(0).nodes().to$().css({cursor: 'pointer'});",
       "  }",
       "};",
       "",
       "// display the child table on click",
       "table.on('click', 'td.details-control', function(){",
       "  var tbl = $(this).closest('table'),",
       "      tblId = tbl.attr('id'),",
       "      td = $(this),",
       "      row = $(tbl).DataTable().row(td.closest('tr')),",
       "      rowIdx = row.index();",
       "  if(row.child.isShown()){",
       "    row.child.hide();",
       "    td.html('&oplus;');",
       "  } else {",
       "    var childId = tblId + '-child-' + rowIdx;",
       "    row.child(format(row.data(), childId)).show();",
       "    td.html('&CircleMinus;');",
       "    format_datatable(row.data(), childId);",
       "  }",
       "});") -> callback

    # return a datatable
    datatable(data=data,
              callback=callback,
              escape=-2, 
              options=list(columnDefs=list(list(visible=FALSE, targets=2:3),
                                           list(orderable=FALSE, className='details-control', targets=1),
                                           list(className='dt-center', targets='_all')),
                           language=list(info='Showing _START_ to _END_ gene modules of _TOTAL_',
                                         infoEmpty='No gene modules to show!'),
                           lengthMenu = list(c(5, 10, 20, -1), 
                                             c(5, 10, 20, 'All'))))
  }) -> output$genes_in_module
}
