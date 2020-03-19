# seuratvis

visualise data contained in a Seurat object

This package is very early stages and under development. Bug reports and feature suggestions are welcome!

# run package

The package expects the current workspace to contain:

* at least one `Seurat` object in the global- or a sub-environment. The global environment and any environments within global are searched for Seurat object which can be selected from the configure tab.

## clone this repository

Make a local copy of the `master` branch using `git clone ChristopherBarrington/seuratvis`

## load the package in `R`

Now use `R` to load the package and run the `shiny` app (you will need to install all of the dependencies manually, sorry!).

```
setwd('/path/to/seuratvis/containing/directory')
devtools::load_all('seuratvis', export_all=FALSE)
launchApp()
```
## install the package

Could use `devtools` to install the master branch as a package. (untested)

```
devtools::install_github('ChristopherBarrington/seuratvis')
launchApp()
```

The app should now open a window using the gene names provided by the vector and the data in the `seurat` object.

# developments

* loading seurat objects by uploading a dataset (maybe)
* all the stuff in the Issues
