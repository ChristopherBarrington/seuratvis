#'
#' 
error_alert <- function(title='Alert title', text='Alert text', session=get(x='session', envir=parent.frame()), btn_labels='OK')
  sendSweetAlert(session=session, title=title, text=text, html=TRUE, closeOnClickOutside=TRUE, type='error')
