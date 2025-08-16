from tempfile import TemporaryDirectory
import time
from typing import Callable
import numpy as np
from scipy.cluster.hierarchy import linkage, fcluster, dendrogram
import pandas as pd
from collections import defaultdict
from scipy.spatial.distance import squareform, pdist
from scipy.cluster import hierarchy
from sklearn.manifold import MDS
import os
import json
import sys
from Bio import SeqIO, Seq, SeqRecord
from workflow.models import RunSettings, WorkflowResult
from transformations import read_csv_matrix
from joblib import Memory


is_testing = "unittest" in sys.modules or any("test" in arg for arg in sys.argv)

if is_testing:
    memory = Memory(location=None, verbose=0)
else:  # pragma: no cover
    memory = Memory(TemporaryDirectory().name, verbose=1)


def run(
    result: WorkflowResult, settings: RunSettings, set_progress: Callable[[int], None]
) -> WorkflowResult:
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

    set_progress(25)

    # TODO: unused - remove?
    dendro = dendrogram(Z, no_plot=True)

    # optimally organize leave indicies for heatmap
    reordered_indices = hierarchy.leaves_list(
        hierarchy.optimal_leaf_ordering(Z, matrix_np)
    )

    set_progress(50)

    # Reorder the matrix
    reordered_matrix = matrix_np[np.ix_(reordered_indices, reordered_indices)]

    set_progress(75)

    if ordered_ids:
        reordered_ids = [ordered_ids[i] for i in reordered_indices]
        reordered_matrix = matrix_np[np.ix_(reordered_indices, reordered_indices)]
    else:
        reordered_matrix = matrix_np[np.ix_(reordered_indices, reordered_indices)]
        reordered_ids = reordered_indices

    set_progress(100)

    return result._replace(
        distance_matrix=reordered_matrix, reordered_ids=reordered_ids
    )


@memory.cache
def calculate_linkage(distance_matrix: np.ndarray, method: str) -> np.ndarray:

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


@memory.cache
def get_mds_coords(distance_mat):
    return MDS(
        n_components=2, dissimilarity="precomputed", random_state=42
    ).fit_transform(distance_mat)


def get_linkage(data: np.ndarray, method: str) -> np.ndarray:
    return calculate_linkage(data, method)


def get_linkage_method_order(data, method, index, threshold=None):
    Z = get_linkage(data, method)

    if threshold is not None:

        cutby = 100 - threshold
        clusters = fcluster(Z, t=cutby, criterion="distance")

        # Create a dictionary mapping cluster ID to sequence indices
        cluster_to_indices = defaultdict(list)
        for i, cluster_id in enumerate(clusters):
            cluster_to_indices[cluster_id].append(i)

        # For each cluster, get the dendrogram order within that cluster
        dendro = dendrogram(Z, no_plot=True)
        leaf_order = dendro["leaves"]

        # create a mapping of the index and its position in the dendrogram
        dendro_position = {idx: pos for pos, idx in enumerate(leaf_order)}

        # sort the clusters by the minimum dendrogram position of their members
        sorted_clusters = sorted(
            cluster_to_indices.keys(),
            key=lambda c: min(dendro_position[idx] for idx in cluster_to_indices[c]),
        )

        # Within each cluster, maintain the dendrogram order
        new_indices = []
        for cluster_id in sorted_clusters:
            cluster_indices = cluster_to_indices[cluster_id]
            # Sort indices within cluster by their dendrogram position
            cluster_indices.sort(key=lambda idx: dendro_position[idx])
            new_indices.extend(cluster_indices)

        # Convert indices to lsit of sequence IDs
        new_order = [index[i] for i in new_indices]
    else:
        #  no threshold,  use dendrogram order
        dendro = dendrogram(Z, no_plot=True)
        leaf_indices = dendro["leaves"]
        new_order = [index[i] for i in leaf_indices]

    return new_order


def get_clusters_dataframe(data, method, threshold, index):
    Z = get_linkage(data, method)
    cluster_data = get_cluster_data_dict(Z, threshold, index)
    return cluster_data_to_dataframe(cluster_data, threshold)


## export fasta files for each clusterr
def export(matrix_path, cluster_data_output_dir, seq_dict_path, threshold, method):
    os.makedirs(cluster_data_output_dir, exist_ok=True)

    ###Prepare cluster.csv output
    cluster_csv_path = os.path.join(cluster_data_output_dir, "cluster.csv")
    with open(matrix_path, "r") as temp_f:
        col_count = [len(l.split(",")) for l in temp_f.readlines()]
        column_names = [i for i in range(0, max(col_count))]

    df = read_csv_matrix(matrix_path)
    index = df.index.tolist()
    data = df.to_numpy()
    data = np.round(data, 2)
    ## save clsuters.csv
    df_result = get_clusters_dataframe(data, method, threshold, index)
    df_result.to_csv(cluster_csv_path, index=False)
    ## Prepare sequence dictionary for fasta export
    seq_dict = {}
    if os.path.exists(seq_dict_path):
        with open(seq_dict_path, "r") as f:
            seq_dict = json.load(f)
    # id the clsuters
    cluster_col_name = [col for col in df_result.columns if "Cluster" in col][0]
    # gather the data
    df_result[cluster_col_name] = df_result[cluster_col_name].astype(int)
    # group it
    grouped_ids_by_cluster = (
        df_result.groupby(cluster_col_name)["ID"].apply(list).to_dict()
    )
    # Create fasta files for each cluster
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


#
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
