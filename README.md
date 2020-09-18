# seuratvis

Visualise data contained in a Seurat object

This package is very early stages and under development. Bug reports and feature suggestions are welcome!

# install package

Either install the package or clone a version for development

## install the package

Could use `remotes` to install the master branch as a package. This method _should_ take care of the dependencies.

```r
remotes::install_github(repo='ChristopherBarrington/seuratvis')
```

Or use the `dev` branch if you are brave.

```r
remotes::install_github(repo='ChristopherBarrington/seuratvis', ref='dev')
```

### possible errors

If a package binary is unavailable for your installed `R` version, you may get an error something like:

```r
Error: (converted from warning) package 'shinydashboardPlus' was built under R version 3.6.3
```

This can be avoided with:

```r
install.packages('shinydashboardPlus', type='source')
```

`dqshiny` may not be available on CRAN so should be installed with the following:

```r
remotes::install_github('daqana/dqshiny')
```

## clone this repository

Make a local copy of the repository using `git clone ChristopherBarrington/seuratvis`. Now use `R` to load the package and run the `shiny` app (you will need to install all of the dependencies manually, sorry!).

```r
setwd('/path/to/seuratvis/containing/directory')
require(devtools) || install.packages('devtools')
devtools::load_all('seuratvis', export_all=FALSE)
```

# run package

The package expects the current workspace to contain:

* at least one `Seurat` object in the global- or a sub-environment. The global environment and any environments within global are searched for Seurat object which can be selected from the configure tab.
  * it will not be happy if there is not at least one!

```r
seuratvis()
```

The app should now open a new window.

# developments

* loading seurat objects by uploading a dataset (maybe)
* all the stuff in the Issues
