#'
#' 
error_alert <- function(session, title, text, btn_labels='OK')
  sendSweetAlert(session=session, title=title, text=text, html=TRUE, closeOnClickOutside=TRUE, type='error')
