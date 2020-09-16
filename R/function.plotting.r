
10^(0:9) -> major_breaks_log10

(2:9) * rep(major_breaks_log10, each=8) -> minor_breaks_log10

#'
#' 
missing_data_plot <- function(label='Data for this plot is unavailable')
  ggplot()+aes()+annotate(geom='text', label=label, x=0, y=0)+theme_void()
