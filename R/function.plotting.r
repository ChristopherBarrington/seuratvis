
10^(0:9) -> major_breaks_log10

(2:9) * rep(major_breaks_log10, each=8) -> minor_breaks_log10

#'
#' 
missing_data_plot <- function(label='Nothing to see here')
  ggplot()+aes()+annotate(geom='text', label=label, x=0, y=0)+theme_void()