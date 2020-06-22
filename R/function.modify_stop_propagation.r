
# https://community.rstudio.com/t/shinydashboard-keep-sidebar-tab-expanded-while-other-tab-is-clicked-expanded/10192/2
modify_stop_propagation <- function(x) {
  x$children[[1]]$attribs$onclick = "event.stopPropagation()"
  x
}
