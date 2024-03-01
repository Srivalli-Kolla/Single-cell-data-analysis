##Running SCVI tools in R## all required packages are installed in conda env scvi-env

#Importing Packages
library(reticulate)
library(anndata)
library(sceasy)
library(Seurat)
library(SeuratData)

#Importing scanpy via py_to_r() function
sc <- import('scanpy', convert = FALSE)  #False will not convert all py outputs to r

#Data loading
data("pbmc3k")
pbmc <- pbmc3k

#Data viewing
pbmc

#Data conversion from seurat to anndata to make it scvi-tools compatible
adata <- convertFormat(pbmc, from="seurat", to="anndata", main_layer="counts", drop_single_values=FALSE)
adata

#To view anndata obs
adata$obs$head()

#To view the python class of data
class(adata$obs)

#Converting obs from py to r
head(py_to_r(adata$obs))

#To view the R class of data
class(py_to_r(adata))

# Convert adata object to R AnnDataR6 object.
adata <- py_to_r(adata)

#Data normalization
X_norm <- sc$pp$normalize_total(adata, target_sum = 1e+09)["X"]
adata$layers["X_norm"] <- X_norm

#Viewing normalized data
head(as.data.frame(adata$layers["X_norm"]))