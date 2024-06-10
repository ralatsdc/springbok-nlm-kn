import cellxgene_census

census = cellxgene_census.open_soma(census_version="latest")

datasets = census["census_info"]["datasets"].read().concat().to_pandas()

lung_obs = (
    census["census_data"]["homo_sapiens"]
    .obs.read(value_filter="tissue_general == 'lung' and is_primary_data == True")
    .concat()
    .to_pandas()
)

lung_datasets = datasets[datasets["dataset_id"].isin(lung_obs["dataset_id"])]
