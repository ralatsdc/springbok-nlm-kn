#!/usr/bin/env python
# coding: utf-8

# # Exploring all data from a tissue
# 
# This tutorial provides a series of examples for how to explore and query the Census in the context of a single tissue, lung. We will summarize cell and gene metadata, then fetch the single-cell expression counts and perform some basic data explorations via [Scanpy](https://scanpy.readthedocs.io/en/stable/) 
# 
# See: https://chanzuckerberg.github.io/cellxgene-census/notebooks/analysis_demo/comp_bio_explore_and_load_lung_data.html
# 
# **Contents**
# 
# 1. Learning about the human lung data.
#    1. Learning about cells of the lung.
#    2. Learning about genes of the lung .
# 2. Fetching all single-cell human lung data from the Census.
# 3. Calculating QC metrics of the lung data.
# 4. Creating a normalized expression layer and embeddings.
# 
# ⚠️ Note that the Census RNA data includes duplicate cells present across multiple datasets. Duplicate cells can be filtered in or out using the cell metadata variable `is_primary_data` which is described in the [Census schema](https://github.com/chanzuckerberg/cellxgene-census/blob/main/docs/cellxgene_census_schema.md#repeated-data).
# 
# ## Learning about the lung data in the Census
# 
# First we will open the Census. If you are not familiar with the basics of the Census API you should take a look at notebook [Learning about the CZ CELLxGENE Census](https://cellxgene-census.readthedocs.io/en/latest/notebooks/analysis_demo/comp_bio_census_info.html)
# 

# In[3]:


import cellxgene_census
import numpy as np
import pandas as pd
import scanpy as sc

census = cellxgene_census.open_soma()


# Let's first take a look at the number of cells from human lung:
# 

# In[4]:


summary_table = census["census_info"]["summary_cell_counts"].read().concat().to_pandas()

summary_table.query("organism == 'Homo sapiens' & category == 'tissue_general' & label =='lung'")


# There you can see the total of cells of under `total_cell_count` and the unique number cells under `unique_cell_count` (i.e. after removing cells that were included in multiple datasets).
# 
# Let's now take a look at the cell and gene information of this slice of the Census.
# 
# ### Learning about cells of lung data
# 
# Let's load the cell metadata for all lung cells and select only the unique cells using `is_primary_data`.
# 

# In[5]:


lung_obs = (
    census["census_data"]["homo_sapiens"]
    .obs.read(value_filter="tissue_general == 'lung' and is_primary_data == True")
    .concat()
    .to_pandas()
)
lung_obs


# You can see that the number or rows represents the total number of unique lung cells in the Census. Now let's take a deeper dive into the characteristics of these cells.
# 
# #### Datasets
# 
# First let's start by looking at what are the datasets and collections from [CELLxGENE Discover](https://cellxgene.cziscience.com/collections) contributing to lung. For this we will use the dataset table at `census["census-info"]["datasets"]` that contains metadata of all datasets used to build this Census.
# 

# In[6]:


census_datasets = (
    census["census_info"]["datasets"]
    .read(column_names=["collection_name", "dataset_title", "dataset_id", "soma_joinid"])
    .concat()
    .to_pandas()
)
census_datasets = census_datasets.set_index("dataset_id")
census_datasets


# The `obs` cell metadata `pandas.DataFrame` contains a column `dataset_id` that can be used for joining to the `census_dataset` `pandas.DataFrame` we just created.
# 
# So let's take a look at the cell counts per `dataset_id` of the lung slice and then join to the dataset table to append the human-readable labels.
# 

# In[7]:


dataset_cell_counts = pd.DataFrame(lung_obs[["dataset_id"]].value_counts())
dataset_cell_counts = dataset_cell_counts.rename(columns={0: "cell_counts"})
dataset_cell_counts = dataset_cell_counts.merge(census_datasets, on="dataset_id")

dataset_cell_counts


# These are all the datasets lung cells whose counts are reprensented in the column `cell_counts`. The top collections with lung data are:
# 
# 1. [The integrated Human Lung Cell Atlas](https://cellxgene.cziscience.com/collections/6f6d381a-7701-4781-935c-db10d30de293).
# 2. [A human cell atlas of fetal gene expression](https://cellxgene.cziscience.com/collections/c114c20f-1ef4-49a5-9c2e-d965787fb90c).
# 3. [High-resolution single-cell atlas reveals diversity and plasticity of tumor-associated neutrophils in non-small cell lung cancer](https://cellxgene.cziscience.com/collections/edb893ee-4066-4128-9aec-5eb2b03f8287).
# 4. [HTAN MSK - Single cell profiling reveals novel tumor and myeloid subpopulations in small cell lung cancer](https://cellxgene.cziscience.com/collections/62e8f058-9c37-48bc-9200-e767f318a8ec).
# 5. [A human fetal lung cell atlas uncovers proximal-distal gradients of differentiation and key regulators of epithelial fates.](https://cellxgene.cziscience.com/collections/2d2e2acd-dade-489f-a2da-6c11aa654028).
# 
# #### Assays
# 
# Let's use similar logic to take a look at all the assays available for human lung data. This tells us that most assays are from 10x technologies and sci-RNA-seq.
# 

# In[8]:


lung_obs[["assay"]].value_counts()


# #### Disease
# 
# And now let's take a look at diseased cell counts, with `normal` indicating non-diseased cells.
# 

# In[9]:


lung_obs[["disease"]].value_counts()


# #### Sex
# 
# There doesn't seem to be strong biases for sex.
# 

# In[10]:


lung_obs[["sex"]].value_counts()


# #### Cell vs nucleus
# 
# The majority of data are from cells and not nucleus.
# 

# In[11]:


lung_obs[["suspension_type"]].value_counts()


# #### Cell types
# 
# Let's take a look at the counts of the top 20 cell types.
# 

# In[12]:


lung_obs[["cell_type"]].value_counts().head(20)


# #### Sub-tissues
# 
# We can look at the original tissue annotations that were mapped to "lung".
# 

# In[13]:


lung_obs[["tissue"]].value_counts()


# ### Learning about genes of lung data
# 
# Let's load the gene metadata of the Census.
# 

# In[14]:


lung_var = census["census_data"]["homo_sapiens"].ms["RNA"].var.read().concat().to_pandas()
lung_var


# You can see the total number of genes represented by the number of rows. This number is actually misleading because it is the join of all genes in the Census. However we know that the lung data comes from a subset of datasets.
# 
# So let's take a look at the number of genes that were measured in each of those datasets.
# 
# To accomplish this we can use the "dataset presence matrix" at `census["census_data"]["homo_sapiens"].ms["RNA"]["feature_dataset_presence_matrix"]`. This is a boolean matrix `N x M` where `N` is the number of datasets and `M` is the number of genes in the Census.
# 
# So we can select the rows corresponding to the lung datasets and perform a row-wise sum.
# 

# In[15]:


presence_matrix = cellxgene_census.get_presence_matrix(census, "Homo sapiens", "RNA")
presence_matrix = presence_matrix[dataset_cell_counts.soma_joinid, :]


# In[16]:


presence_matrix.sum(axis=1).A1


# In[17]:


genes_measured = presence_matrix.sum(axis=1).A1
dataset_cell_counts["genes_measured"] = genes_measured
dataset_cell_counts


# You can see the genes measured in each dataset represented in `genes_measured`. Now lets get the **genes that were measured in all datasets**.
# 

# In[18]:


var_somaid = np.nonzero(presence_matrix.sum(axis=0).A1 == presence_matrix.shape[0])[0].tolist()


# In[19]:


lung_var = lung_var.query(f"soma_joinid in {var_somaid}")
lung_var


# The number of rows represents the genes that were measured in all lung datasets.
# 
# ### Summary of lung metadata
# 
# In the previous sections, using the Census we learned the following information:
# 
# - The total number of unique lung cells and their composition for:
#   - Number of datasets.
#   - Number sequencing technologies, most of which are 10x
#   - Mostly human data, but some diseases exist, primarily "lung adenocarcinoma" and "COVID-19 infected"
#   - No sex biases.
#   - Mostly data from cells (\~80%) rather than nucleus (\~20%)
# - A total of **~12k** genes were measured across all cells.
# 
# ##  Fetching all single-cell human lung data from the Census
# 
# Since loading the entire lung data is resource-intensive, for the sake of this exercise let's load a subset of the lung data into an `anndata.AnnData` object and perform some exploratory analysis. 
# 
# We will subset to 100,000 random unique cells using the `lung_obs` `pandas.DataFrame` we previously created.

# In[20]:


lung_cell_subsampled_n = 100
lung_cell_subsampled_ids = lung_obs["soma_joinid"].sample(lung_cell_subsampled_n, random_state=1).tolist()


# Now we can directly use the values of `soma_joinid` for querying the Census data and obtaining an `AnnData` object.

# In[21]:


lung_gene_ids = lung_var["soma_joinid"].to_numpy()
lung_adata = cellxgene_census.get_anndata(
    census,
    organism="Homo sapiens",
    obs_coords=lung_cell_subsampled_ids,
    var_coords=lung_gene_ids,
)

lung_adata.var_names = lung_adata.var["feature_name"]


# In[24]:


lung_adata


# We are done with the census, so close it

# In[25]:


census.close()
del census


# lung_data_saved = lung_data.copy()

# ## Calculating QC metrics of the lung data
# 
# Now let's take a look at some QC metrics
# 
# **Top genes per cell**
# 

# In[22]:


sc.pl.highest_expr_genes(lung_adata, n_top=20)


# **Number of sequenced genes by assay**
# 

# In[23]:


sc.pp.calculate_qc_metrics(lung_adata, percent_top=None, log1p=False, inplace=True)
sc.pl.violin(lung_adata, "n_genes_by_counts", groupby="assay", jitter=0.4, rotation=90)


# **Total counts by assay**
# 

# In[24]:


sc.pl.violin(lung_adata, "total_counts", groupby="assay", jitter=0.4, rotation=90)


# You can see that Smart-Seq2 is an outlier for the total counts per cell, so let's exlcude it to see how the rest of the assays look like
# 

# In[25]:


sc.pl.violin(
    lung_adata[lung_adata.obs["assay"] != "Smart-seq2",],
    "total_counts",
    groupby="assay",
    jitter=0.4,
    rotation=90,
)


# ## Creating a normalized expression layer and embeddings
# 
# Let's perform a bread and butter normalization and take a look at UMAP embeddings, but for all the data below we'll exclude Smart-seq2 as this requires an extra step to normalize based on gene lengths
# 

# In[26]:


lung_adata = lung_adata[lung_adata.obs["assay"] != "Smart-seq2",].copy()
lung_adata.layers["counts"] = lung_adata.X


# Now let's do some basic normalization:
# 
# - Normalize by sequencing depth
# - Transform to log-scale
# - Select 500 highly variable genes
# - Scale values across the gene axis
# 

# In[27]:


sc.pp.normalize_total(lung_adata, target_sum=1e4)
sc.pp.log1p(lung_adata)
sc.pp.highly_variable_genes(lung_adata, n_top_genes=500, flavor="seurat_v3", layer="counts")
lung_adata = lung_adata[:, lung_adata.var.highly_variable]
sc.pp.scale(lung_adata, max_value=10)


# And reduce dimensionality by obtaining UMAP embeddings.
# 

# In[28]:


sc.tl.pca(lung_adata)
sc.pp.neighbors(lung_adata)
sc.tl.umap(lung_adata)


# And plot these embeddings.
# 

# In[29]:


n_cell_types = len(lung_adata.obs["cell_type"].drop_duplicates())

from random import randint

colors = []

for i in range(len(lung_adata.obs["cell_type"].drop_duplicates())):
    colors.append("#%06X" % randint(0, 0xFFFFFF))


# In[30]:


sc.pl.umap(lung_adata, color="cell_type", palette=colors, legend_loc=None)


# Let's color by assay.
# 

# In[31]:


sc.pl.umap(lung_adata, color="assay")


# Given the high number of cell types it makes it hard to visualize, so let's look at the top 20 most abundant cell types.
# 

# In[32]:


top_cell_types = lung_adata.obs["cell_type"].value_counts()
top_cell_types = list(top_cell_types.reset_index().head(20)["cell_type"])


# In[33]:


lung_adata_top_cell_types = lung_adata[[i in top_cell_types for i in lung_adata.obs["cell_type"]], :]
sc.pl.umap(lung_adata_top_cell_types, color="cell_type")


# Let's color by assay of this subset of the data.
# 

# In[34]:


sc.pl.umap(lung_adata_top_cell_types, color="assay")


# In[ ]:





# In[ ]:




