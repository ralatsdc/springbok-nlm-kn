import os

import nsforest as ns
from nsforest import nsforesting
import scanpy as sc

DATA_DIR = "../data"

CELLXGENE_DIR = f"{DATA_DIR}/cellxgene"

NSFOREST_DIR = f"{DATA_DIR}/nsforest-2024-06-27"
TOTAL_COUNTS = 5000  # TODO: Select a more sensible value

NCBI_CELL_DIR = f"{DATA_DIR}/ncbi-cell"


def run_nsforest_on_file(h5ad_filename, cluster_header="cell_type"):
    """Run NSForest using the specified dataset filename, and
    cluster_header.

    Parameters
    ----------
    h5ad_filename : str
       The dataset filename
    cluster_header : str
       The cluster header

    Returns
    -------
    None
    """
    # Assign results filename and directory
    pp_h5ad_filename = f"pp_{h5ad_filename}"
    results_dirname = h5ad_filename.split(".")[0]
    results_dirpath = f"{NSFOREST_DIR}/{results_dirname}"

    # Run NSForest if results do not exist
    if not os.path.exists(results_dirpath):
        os.makedirs(results_dirpath)

        print(f"Loading unprocessed AnnData file: {h5ad_filename}")
        h5ad_filepath = f"{CELLXGENE_DIR}/{h5ad_filename}"
        up_adata = sc.read_h5ad(h5ad_filepath)

        # TODO: Check validity of downsampling
        print("Calculating QC metrics")
        up_metrics = sc.pp.calculate_qc_metrics(up_adata)
        if up_metrics[1]["total_counts"].sum() > TOTAL_COUNTS:
            print("Downsampling unprocessed AnnData file")
            ds_adata = sc.pp.downsample_counts(
                up_adata, total_counts=TOTAL_COUNTS, copy=True
            )
        else:
            ds_adata = up_adata  # No need to copy

        print("Generating scanpy dendrogram")
        # Dendrogram order is stored in
        # `pp_adata.uns["dendrogram_cluster"]["categories_ordered"]`
        pp_adata = up_adata.copy()
        pp_adata.obs[cluster_header] = pp_adata.obs[cluster_header].astype(str)
        pp_adata.obs[cluster_header] = pp_adata.obs[cluster_header].astype("category")
        pp_adata = ns.pp.dendrogram(
            pp_adata,
            cluster_header,
            save=False,
            output_folder=results_dirpath,
            outputfilename_suffix=cluster_header,
        )

        print("Calculating cluster medians per gene")
        pp_adata = ns.pp.prep_medians(pp_adata, cluster_header)

        print("Calculating binary scores per gene per cluster")
        pp_adata = ns.pp.prep_binary_scores(pp_adata, cluster_header)

        pp_h5ad_filepath = f"{results_dirpath}/{pp_h5ad_filename}"
        print(f"Saving preprocessed AnnData file: {pp_h5ad_filepath}")
        pp_adata.write_h5ad(pp_h5ad_filepath)

        print(f"Running NSForest for preprocessed AnnData file: {pp_h5ad_filename}")
        results = nsforesting.NSForest(
            pp_adata,
            cluster_header,
            output_folder=f"{results_dirpath}/",
            outputfilename_prefix=cluster_header,
        )

    else:
        print(f"Completed NSForest for preprocessed AnnData file: {pp_h5ad_filename}")
