import cellxgene_census
import numpy as np
import pandas as pd
import scanpy as sc

census = cellxgene_census.open_soma()


# A data frame with high-level information of this Census, e.g. build date, total cell count, etc.
# Index(['soma_joinid', 'label', 'value'], dtype='object')
summary = census["census_info"]["summary"].read().concat().to_pandas()

# A data frame with all datasets from CELLxGENE Discover used to create the Census.
# Index(['soma_joinid', 'citation', 'collection_id', 'collection_name',
#        'collection_doi', 'dataset_version_id', 'dataset_title',
#        'dataset_h5ad_path', 'dataset_total_cell_count'],
#       dtype='object')
datasets = census["census_info"]["datasets"].read().concat().to_pandas()
datasets = datasets.set_index("dataset_id")

# A data frame with cell counts stratified by relevant cell metadata
# Index(['soma_joinid', 'organism', 'category', 'label', 'ontology_term_id',
#        'total_cell_count', 'unique_cell_count'],
#       dtype='object')
counts = census["census_info"]["summary_cell_counts"].read().concat().to_pandas()

# Data matrices, currently only raw counts exist X["raw"]
# X = census["census_data"]["homo_sapiens"].ms["RNA"].X
adata = cellxgene_census.get_anndata(
    census=census,
    organism="Homo sapiens",
)

# Cell metadata
# obs = census["census_data"]["homo_sapiens"].obs
# Index(['soma_joinid', 'dataset_id', 'assay', 'assay_ontology_term_id',
#        'cell_type', 'cell_type_ontology_term_id', 'development_stage',
#        'development_stage_ontology_term_id', 'disease',
#        'disease_ontology_term_id', 'donor_id', 'is_primary_data',
#        'observation_joinid', 'self_reported_ethnicity',
#        'self_reported_ethnicity_ontology_term_id', 'sex',
#        'sex_ontology_term_id', 'suspension_type', 'tissue',
#        'tissue_ontology_term_id', 'tissue_type', 'tissue_general',
#        'tissue_general_ontology_term_id', 'raw_sum', 'nnz', 'raw_mean_nnz',
#        'raw_variance_nnz', 'n_measured_vars'],
#       dtype='object')
obs = cellxgene_census.get_obs(census, "homo_sapiens")

# Gene Metadata
# var = census["census_data"]["homo_sapiens"].ms["RNA"].var
# Index(['soma_joinid', 'feature_id', 'feature_name', 'feature_length', 'nnz',
#        'n_measured_obs'],
#       dtype='object')
var = cellxgene_census.get_var(census, "homo_sapiens")

# See: https://chanzuckerberg.github.io/cellxgene-census/notebooks/api_demo/census_query_extract.html

# Convenience wrapper around tiledbsoma.Experiment query, to build and execute a query, and return it as an anndata.AnnData object
# https://chanzuckerberg.github.io/cellxgene-census/_autosummary/cellxgene_census.get_anndata.html#cellxgene_census.get_anndata
# cellxgene_census.get_anndata
# adata = cellxgene_census.get_anndata(
#     census=census,
#     organism="Homo sapiens",
#     var_value_filter="feature_id in ['ENSG00000161798', 'ENSG00000188229']",
#     obs_value_filter="cell_type == 'B cell' and tissue_general == 'lung' and disease == 'COVID-19' and is_primary_data == True",
#     obs_column_names=["sex"],
# )

# Get the observation metadata for a query on the census.
# https://chanzuckerberg.github.io/cellxgene-census/_autosummary/cellxgene_census.get_obs.html#cellxgene_census.get_obs
# cellxgene_census.get_obs
# lung_obs = cellxgene_census.get_obs(
#     census, "homo_sapiens", value_filter="tissue_general == 'lung' and is_primary_data == True"
# )

# Get the variable metadata for a query on the census.
# https://chanzuckerberg.github.io/cellxgene-census/_autosummary/cellxgene_census.get_var.html#cellxgene_census.get_var
# cellxgene_census.get_var
# lung_var = cellxgene_census.get_var(
#     census, "homo_sapiens", value_filter="tissue_general == 'lung' and is_primary_data == True"
# )

census.close()