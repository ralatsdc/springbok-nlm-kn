
process preprocess_cxg_file {

    input:
    path cxg_file
    val cluster_header
    val results_header

    output:
    path "*_pp.h5ad"
    path "clusters.txt"
    path "collected_${cluster_header}_results.csv"
    
    script:
    """
    # Preprocess the specified file
    nsforest.py --preprocess-adata-file -c ${cluster_header} ${cxg_file}

    # Write header to results file
    echo ${results_header} > "collected_${cluster_header}_results.csv"
    """

}

process run_nsforest_by_cluster {

    input:
    tuple path(pp_cxg_file), val(cluster)
    val cluster_header

    output:
    path "${cluster_header}_${cluster}_results.csv"
    
    script:
    """
    # Run NSforest
    nsforest.py --run-nsforest-without-preprocessing -c ${cluster_header} -l ${cluster} ${pp_cxg_file}
    cp "${cluster_header}_results.csv" "${cluster_header}_${cluster}_results.csv"
    """
}

process collect_results {

    input:
    path results_file
    val results_header
    path collected_results_file

    output:
    path collected_results_file

    script:
    """
    cat ${results_file} | grep -v ${results_header} >> ${collected_results_file}
    """

}

workflow {

    cxg_files_ch = channel.fromPath("data/*.h5ad")
    cluster_header_ch = channel.value("cluster")
    results_header_ch = channel.value(
        "clusterName,clusterSize,f_score,PPV,recall,TN,FP,FN,TP,marker_count,NSForest_markers,binary_genes,onTarget"
    )

    (
        pp_cxg_files_ch, cluster_names_file_ch, collected_results_file_ch
    ) = preprocess_cxg_file(
        cxg_files_ch, cluster_header_ch, results_header_ch
    )

    pp_cxg_file_and_clusters_ch = pp_cxg_files_ch.combine(
        cluster_names_file_ch.splitText(){ it.trim() }
    )

    results_files_ch = run_nsforest_by_cluster(
        pp_cxg_file_and_clusters_ch, cluster_header_ch
    )

    final_results_file_ch = collect_results(
        results_files_ch.toList(), results_header_ch, collected_results_file_ch
    )

    final_results_file_ch.view()
}
