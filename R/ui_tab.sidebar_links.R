
.ui_tab.sidebar_links <- function() {
  bquote({
    # menu tab hyperlinks
    list(hr=hr(),
         email_me=menuItem(text='mail me', href='mailto:christopher.barrington@crick.ac.uk?subject=[seuratvis] Hello there', icon=icon('comment-dots'), newtab=FALSE),
         bug_report=menuItem(text='report a bug', href='mailto:christopher.barrington@crick.ac.uk?subject=[seuratvis] I found a bug!', icon=icon('bug'), newtab=FALSE),
         feature_request=menuItem(text='suggest a feature', href='mailto:christopher.barrington@crick.ac.uk?subject=[seuratvis] It would be cool if...', icon=icon('lightbulb'), newtab=FALSE),
         github_link=menuItem(text=sprintf('GitHub (version: %s)', packageVersion('seuratvis')), href='github.com/ChristopherBarrington/seuratvis', icon=icon('code-branch'), newtab=FALSE)) %>%
      append(x=menus) -> menus})
}
