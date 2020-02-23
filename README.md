# seuratvis

visualise data contained in a Seurat object

This package is very early stages and under development. Bug reports and feature suggestions are welcome!

# run package

At the moment this is not very strightforward and needs to be (dramatically) improved. The package expects the current workspace to contain:

* a `Seurat` object named `seurat`
* (optionally) a named `character` vector called `gene_names_to_description`, where names are gene IDs (ie the row names of `seurat`) and the values are a description of the gene

## clone this repository

Make a local copy of the `master` branch using `git clone ChristopherBarrington/seuratvis`

## load the package in `R`

Now use `R` to load the package and run the `shiny` app

```
setwd('/path/to/seuratvis/containing/directory')
devtools::load_all('seuratvis', export_all=FALSE)
launchApp()
```

The app should now open a window using the gene names provided by the vector and the data in the `seurat` object.

# developments

* installation process - this will become an installation using `devtools::install_github` instead.
* loading seurat objects by uploading a dataset (maybe)
* check for `seurat` object before continuing installation etc or remove dependence for the object to be in workspace for installation
