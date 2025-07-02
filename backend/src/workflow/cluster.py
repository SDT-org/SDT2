from tempfile import TemporaryDirectory
import time
import numpy as np
from scipy.cluster.hierarchy import linkage, fcluster, dendrogram
import pandas as pd
from collections import defaultdict
from scipy.spatial.distance import squareform, pdist
from scipy.cluster import hierarchy
from sklearn.manifold import MDS
from joblib import Memory
import os
import json
from Bio import SeqIO, Seq, SeqRecord
from workflow.models import RunSettings, WorkflowResult
from transformations import read_csv_matrix


def run(result: WorkflowResult, settings: RunSettings) -> WorkflowResult:
    distance_matrix = result.distance_matrix
    if distance_matrix is None:
        return result

    cluster_method = settings.cluster_method
    ordered_ids = result.ordered_ids or []
    # Convert to numpy array if needed
    if isinstance(distance_matrix, np.ndarray):
        matrix_np = distance_matrix
    elif hasattr(distance_matrix, "to_numpy"):
        matrix_np = distance_matrix.to_numpy()
    else:
        matrix_np = np.array(distance_matrix)

    Z = calculate_linkage(matrix_np, cluster_method)

    dendro = dendrogram(Z, no_plot=True)
    ## squareform for optimpal leaf ordering
    Z = hierarchy.ward(matrix_np)
    
    #reordered_indices = dendro["leaves"]
    #organize leave indicies for heatmap
    reordered_indices = hierarchy.leaves_list(hierarchy.optimal_leaf_ordering(Z, matrix_np))
    # Reorder the matrix
    reordered_matrix = matrix_np[np.ix_(reordered_indices, reordered_indices)]

    if ordered_ids:
        reordered_ids = [ordered_ids[i] for i in reordered_indices]
        if hasattr(distance_matrix, "loc"):
            reordered_matrix = matrix_np[np.ix_(reordered_indices, reordered_indices)]
        else:
            reordered_matrix = matrix_np[np.ix_(reordered_indices, reordered_indices)]
    else:
        reordered_matrix = matrix_np[np.ix_(reordered_indices, reordered_indices)]
        reordered_ids = reordered_indices

    return result._replace(
        distance_matrix=reordered_matrix, reordered_ids=reordered_ids
    )


def calculate_linkage(distance_matrix: np.ndarray, method: str) -> np.ndarray:
    """Calculate hierarchical clustering linkage matrix."""
    if method in ["ward", "centroid", "median"]:
        # MDS for methods that need Euclidean distance
        start = time.perf_counter()
        mds_coords = get_mds_coords(distance_matrix)
        if mds_coords is None:
            raise ValueError("MDS failed to compute coordinates.")
        end = time.perf_counter()
        print(f"MDS coords computed in {end - start:.2f} seconds")
        y = pdist(mds_coords)
        metric = "euclidean"
    else:
        # For other methods, use the distance matrix directly
        start = time.perf_counter()
        y = squareform(distance_matrix)
        metric = "precomputed"
        end = time.perf_counter()
        print(f"Distance matrix prepared in {end - start:.2f} seconds")
    # print(method, metric, y)
    return linkage(y, method=method, metric=metric)


def get_mds_coords(distance_mat):
    return MDS(
        n_components=2, dissimilarity="precomputed", random_state=42
    ).fit_transform(distance_mat)


def get_linkage(data: np.ndarray, method: str) -> np.ndarray:
    return calculate_linkage(data, method)


def get_linkage_method_order(data, method, index):
    Z = get_linkage(data, method)

    # Create dendrogram to get leaf order
    dendro = dendrogram(Z, no_plot=True)
    leaf_indices = dendro["leaves"]
    new_order = [index[i] for i in leaf_indices]

    return new_order


def get_clusters_dataframe(data, method, threshold, index):
    Z = get_linkage(data, method)
    cluster_data = get_cluster_data_dict(Z, threshold, index)
    return cluster_data_to_dataframe(cluster_data, threshold)


def export(matrix_path, cluster_data_output_dir, seq_dict_path, threshold, method):
    os.makedirs(cluster_data_output_dir, exist_ok=True)
    cluster_csv_path = os.path.join(cluster_data_output_dir, "cluster.csv")

    with open(matrix_path, "r") as temp_f:
        col_count = [len(l.split(",")) for l in temp_f.readlines()]
        column_names = [i for i in range(0, max(col_count))]

    df = read_csv_matrix(matrix_path)
    index = df.index.tolist()
    data = df.to_numpy()
    data = np.round(data, 2)

    df_result = get_clusters_dataframe(data, method, threshold, index)
    df_result.to_csv(cluster_csv_path, index=False)

    seq_dict = {}
    if os.path.exists(seq_dict_path):
        with open(seq_dict_path, "r") as f:
            seq_dict = json.load(f)

    cluster_col_name = [col for col in df_result.columns if "Cluster" in col][0]
    df_result[cluster_col_name] = df_result[cluster_col_name].astype(int)
    grouped_ids_by_cluster = (
        df_result.groupby(cluster_col_name)["ID"].apply(list).to_dict()
    )
    for cluster_num, list_of_seq_ids in grouped_ids_by_cluster.items():
        records_for_this_cluster = []
        for seq_id in list_of_seq_ids:
            if seq_id in seq_dict:
                sequence_string = seq_dict[seq_id]
                bio_seq = Seq.Seq(sequence_string)
                seq_record_obj = SeqRecord.SeqRecord(
                    bio_seq, id=str(seq_id), description=""
                )
                records_for_this_cluster.append(seq_record_obj)
        if records_for_this_cluster:
            fasta_filename = os.path.join(
                cluster_data_output_dir, f"cluster_{cluster_num}.fasta"
            )
            with open(fasta_filename, "w") as output_handle:
                SeqIO.write(records_for_this_cluster, output_handle, "fasta")

    return df_result


def get_cluster_data_dict(Z, threshold, index):
    # Set cut threshold
    cutby = 100 - threshold
    # Identify clusters from threshold cut
    clusters = fcluster(Z, t=cutby, criterion="distance")
    # Store as dict
    cluster_dict = defaultdict(list)
    for i, label in enumerate(clusters):
        cluster_dict[label - 1].append(index[i])

    return cluster_dict


def cluster_data_to_dataframe(cluster_data: dict, threshold: int):
    flattened_output = [
        (item, key + 1) for key, sublist in cluster_data.items() for item in sublist
    ]
    dataframe = pd.DataFrame(flattened_output)
    dataframe.columns = ["ID", "Cluster - Threshold: " + str(threshold)]

    return dataframe


def order_clusters_sequentially(clusters):
    # Let's make it NP
    labels = np.array(clusters)
    # Create newly reordered labels that make sense
    new_order_labels = np.searchsorted(np.unique(labels), labels)
    return (new_order_labels + 1).tolist()
