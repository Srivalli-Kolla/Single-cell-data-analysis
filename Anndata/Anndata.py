# %% [markdown]
# Annotating data which is in matrix (N x D) 
# **obs** N observations - Samples or cell types - uni dimensional
# **var** D dimensional vectors -  Variables - uni dimensional
# **Layers**  different forms of our original core data - may be normalized data in log form and raw data
# 
#         var var var
#         gene gene gene
# obs cell
# obs cell
# obs cell
# 
# **Metadata**
# 
# 1. Observation/Variable level - single set of information is added on both obs and var side
# 2. Observation/variable-level matrices - More than one set of information is addedd. Usually it will be in form of matrix.
# 
#         varm varm varm
# obsm
# obsm
# obsm
# 
# Both **obsm and varm **are additional obs and var of metadata but not obs and var original h5ad data and are multidimensional
# 
# ## DIMENSIONS OF H5AD DATA AND METADATA MATRIX HAS TO BE SAME ##
# 
# 3. Unstructured metadata - This can be anything, like a list or a dictionary with some general information that was useful in the analysis of our data.
# 
# Square matrices representing graphs are stored in **obsp and varp**, with both of their own dimensions aligned to their associated axis
# 
# obs - One-dimensional annotation of observations (pd.DataFrame).
# obs_names - Names of observations (alias for .obs.index).
# obsm - Multi-dimensional annotation of observations (mutable structured ndarray).
# obsp - Pairwise annotation of observations, a mutable mapping with array-like values.
# 
# var - One-dimensional annotation of variables/ features (pd.DataFrame).
# var_names - Names of variables (alias for .var.index).
# varm - Multi-dimensional annotation of variables/features (mutable structured ndarray).
# varp - Pairwise annotation of variables/features, a mutable mapping with array-like values.

# %%
#importing packages # Make sure that all required packages are downloaded in given environment (conda activate anndata)
import numpy as np
import pandas as pd
import anndata as ad
from scipy.sparse import csr_matrix
print(ad.__version__)

# %%
#INITIALIZING ANNDATA

#Reading h5ad data
h5data = ad.read_h5ad("/Users/srivalli/Desktop/Heart/hca_heart_immune_download.h5ad")

#Data structure of anndata file i.e., summary stastics of the data
ad.AnnData(h5data)
adata = ad.AnnData(h5data)

#making matrix with the anndata
matrix = adata.X


# %%
#SUBSETTING DATA

#Subsetting data when gender is female
bdata  = adata[adata.obs.cell_source == "Female"]
bdata

#To view it
bdata.obs

#Subsetting data when cell_sopu is female
cdata  = adata[adata.obs.scNym == "CD4+T_cell"]
cdata

#To view it
cdata.obs

# %%
#ADDING ALIGNED METADATA   

#Adding more information to the dataset which is aligned 
ct = np.random.choice(["B", "T", "Monocyte"], size=(adata.n_obs,))
adata.obs["cell_type"] = pd.Categorical(ct)  # Categoricals are preferred for efficiency
adata.obs

# %%
#MAKING LAYERS OF DATA

#Normalized data using log transformation
adata.layers["log_transformed"] = np.log1p(adata.X)
adata

#Making a dataframe of log values
logdata = adata.to_df(layer="log_transformed")
logdata

# %%
#Saving data into file
adata.write('my_results.h5ad', compression="gzip")
logdata.to_csv('logdata.txt', compression="gzip")


