import os
import sys
import json
import numpy as np
from pandas import DataFrame
from numpy import eye, where, nan, nanmin, nanmax

sys.path.append(os.path.join(os.path.dirname(__file__), ".."))

from document_paths import build_document_paths
from transformations import read_csv_matrix, read_stats_csv, to_triangle
from file_utils import read_json_file
from workflow import cluster


def file_exists(path):
    return os.path.exists(path)


def get_data(doc_id: str, get_document, update_document):
    doc = get_document(doc_id)
    doc_paths = build_document_paths(doc.tempdir_path)

    df = read_csv_matrix(doc_paths.matrix)
    ids = df.index.tolist()
    update_document(doc_id, sequences_count=len(ids))
    data = df.to_numpy()

    if os.path.exists(doc_paths.stats):
        stats_df = read_stats_csv(doc_paths.stats)
    else:
        stats_df = DataFrame([])

    # Generate identity scores from matrix for CSV imports
    identity_scores = []
    all_scores = []
    id_map = {id: idx for idx, id in enumerate(ids)}
    
    # Check if LZANI was used
    is_lzani = False
    if file_exists(doc_paths.run_settings):
        run_settings = read_json_file(doc_paths.run_settings)
        is_lzani = run_settings.get("analysis_method") == "lzani"
    
    # Matrix contains distance values, convert to similarity for identity scores
    for i in range(len(ids)):
        for j in range(i):
            identity_score = 100 - data[i, j]
            all_scores.append(identity_score)
            
            if not (is_lzani and identity_score < 1):
                identity_scores.append([id_map[ids[i]], id_map[ids[j]], identity_score])
    
    # Calculate unaligned count after collecting all scores
    unaligned_count = len([s for s in all_scores if s < 3]) if is_lzani else 0
    
    # Convert matrix to triangle format for heatmap
    heat_data = DataFrame(data, index=ids)
    heat_data = to_triangle(heat_data)

    # Calculate min/max after conversion to similarity values
    heat_data_np = heat_data.to_numpy()
    diag_mask = eye(heat_data_np.shape[0], dtype=bool)
    heat_data_no_diag = where(diag_mask, nan, heat_data_np)
    
    # For LZ-ANI, adjust minimum to 70 (or lowest value above 70)
    if is_lzani:
        valid_values = heat_data_no_diag[~np.isnan(heat_data_no_diag) & (heat_data_no_diag >= 70)]
        if valid_values.size > 0:
            min_val = max(70, int(nanmin(valid_values)))
        else:
            min_val = 30
    else:
        min_val = int(nanmin(heat_data_no_diag))
    
    max_val = int(nanmax(heat_data_no_diag))

    parsedData = heat_data.values.tolist()

    metadata = dict(minVal=min_val, maxVal=max_val)
    if file_exists(doc_paths.run_settings):
        metadata["run"] = read_json_file(doc_paths.run_settings)
    
    if is_lzani:
        metadata["unaligned_count"] = unaligned_count
    
    # Helper function to calculate statistics
    def calculate_stats(values):
        if len(values) == 0:
            return None
        return {
            "mean": float(np.mean(values)),
            "median": float(np.median(values)),
            "std": float(np.std(values)),
            "min": float(np.min(values)),
            "max": float(np.max(values)),
            "q1": float(np.percentile(values, 25)),
            "q3": float(np.percentile(values, 75)),
            "count": len(values)
        }
    
    if identity_scores:
        scores_array = np.array([score[2] for score in identity_scores])
        metadata["distribution_stats"] = calculate_stats(scores_array)
    
    if not stats_df.empty:
        stats_array = np.array(stats_df.values)
        metadata["gc_stats"] = calculate_stats(stats_array[:, 1])
        metadata["length_stats"] = calculate_stats(stats_array[:, 2])

    data_to_dump = dict(
        metadata=metadata,
        data=[ids] + parsedData,
        ids=ids,
        identity_scores=identity_scores,
        full_stats=stats_df.values.tolist(),
    )
    return json.dumps(data_to_dump)


def get_clustermap_data(doc_id: str, threshold: float, method: str, get_document):
    doc = get_document(doc_id)
    doc_paths = build_document_paths(doc.tempdir_path)
    
    matrix_df = read_csv_matrix(doc_paths.matrix)
    matrix_np = matrix_df.to_numpy()
    matrix_np = np.round(matrix_np, 2)
    sorted_ids = matrix_df.index.tolist()

    # Get cluster assignments
    seqid_clusters_df = cluster.get_clusters_dataframe(
        matrix_np, method, threshold, sorted_ids
    )
    seqid_clusters_df = seqid_clusters_df.rename(
        columns={
            str(seqid_clusters_df.columns[0]): "id",
            str(seqid_clusters_df.columns[1]): "cluster",
        }
    )

    seqid_clusters_df["id"] = seqid_clusters_df["id"].astype(str)
    seqid_clusters_df["original_cluster"] = seqid_clusters_df["cluster"]

    # Get the linkage-based order to group sequences by cluster
    new_order = cluster.get_linkage_method_order(matrix_np, method, sorted_ids, threshold)

    # Create a mapping from old index to new index
    reorder_indices = [sorted_ids.index(id) for id in new_order]

    # Reorder the matrix
    reordered_matrix = matrix_np[np.ix_(reorder_indices, reorder_indices)]

    # Reorder the tick text
    reordered_tick_text = [str(sorted_ids[i]) for i in reorder_indices]

    # Now assign sequential cluster numbers based on visual order
    id_to_cluster = {row["id"]: row["cluster"] for _, row in seqid_clusters_df.iterrows()}

    # get cluster stats
    seen_clusters = {}
    next_cluster_num = 1
    
    total_clusters= len(seqid_clusters_df["cluster"].unique())

    largest_cluster = max(seqid_clusters_df["cluster"].value_counts())
 
    cluster_counts = seqid_clusters_df["cluster"].value_counts()
    singleton_clusters = cluster_counts[cluster_counts == 1]
    singletons = len(singleton_clusters)
    print(singletons, "singletons")
    
    # Go through sequences in descending order from topmost cluster and assign ascending sequential cluster numbers
    for seq_id in reordered_tick_text:
        original_cluster = id_to_cluster.get(seq_id)
        if original_cluster is not None and original_cluster not in seen_clusters:
            seen_clusters[original_cluster] = next_cluster_num
            next_cluster_num += 1

    # Update the dataframe with new sequential cluster numbers
    for idx, row in seqid_clusters_df.iterrows():
        original_cluster = row["cluster"]
        if original_cluster in seen_clusters:
            seqid_clusters_df.at[idx, "cluster"] = seen_clusters[original_cluster]

    # Update the cluster data to match the new order so that the legend matches
    cluster_data = seqid_clusters_df.to_dict(orient="records")

    cluster_stats = {
        "total_clusters": total_clusters,
        "largest_cluster": int(largest_cluster),
        "singleton_clusters": singletons,
    }

    return {
        "matrix": to_triangle(reordered_matrix, fill_value=None).tolist(),
        "tickText": reordered_tick_text,
        "clusterData": cluster_data,
        "cluster_stats": cluster_stats,
    }
