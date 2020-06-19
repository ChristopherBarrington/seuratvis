reformat_project_name <- function(x)
  x %>% str_replace_all(pattern='_', replacement=' ') %>% str_to_upper()

