# seuratvis

Visualise data contained in a Seurat object

This package is very early stages and under development. Bug reports and feature suggestions are welcome!

# install package

Either install the package or clone a version for development

## install the package

Could use `devtools` to install the master branch as a package. This method _should_ take care of the dependencies.

```r
require(devtools) || install.packages('devtools')
install_github('ChristopherBarrington/seuratvis')
```

## clone this repository

Make a local copy of the repository using `git clone ChristopherBarrington/seuratvis`. Now use `R` to load the package and run the `shiny` app (you will need to install all of the dependencies manually, sorry!).

```r
setwd('/path/to/seuratvis/containing/directory')
devtools::load_all('seuratvis', export_all=FALSE)
```

# run package

The package expects the current workspace to contain:

* at least one `Seurat` object in the global- or a sub-environment. The global environment and any environments within global are searched for Seurat object which can be selected from the configure tab.
  * it will not be happy if there is not at least one!

```r
library(seuratvis)
launchApp()
```

The app should now open a new window.

# developments

* loading seurat objects by uploading a dataset (maybe)
* all the stuff in the Issues
