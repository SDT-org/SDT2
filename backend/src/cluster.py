import os
from tempfile import TemporaryDirectory
import numpy as np
from scipy.cluster.hierarchy import linkage, fcluster, dendrogram
import pandas as pd
from collections import defaultdict
from scipy.spatial.distance import squareform, pdist
from sklearn.manifold import MDS
import time
from joblib import Memory
cache_dir = TemporaryDirectory()
memory = Memory(cache_dir.name, verbose=1)

@memory.cache
def calculate_linkage(data, method) -> np.ndarray:
    dist_data = data.copy()  # TODO: do we need to copy the data?

    # Check if it's a similarity matrix (values near 100) or distance matrix
    if np.max(dist_data) > 90:
        distance_mat = 100 - dist_data
    else:
        distance_mat = dist_data

    # MDS for methods that need Euclidean distance
    mds_coords = MDS(
        n_components=2, dissimilarity="precomputed", random_state=42
    ).fit_transform(distance_mat)

    if method in ["ward", "centroid", "median"]:
        y = pdist(mds_coords)
        metric = "euclidean"
    else:
        y = squareform(distance_mat)
        metric = "precomputed"

    return linkage(y, method=method, metric=metric)


def get_linkage(data: np.ndarray, method: str) -> np.ndarray:
    start = time.perf_counter()
    result = calculate_linkage(data, method)
    end = time.perf_counter()
    print(f"Calculate linkage took {end - start:0.4f} seconds")
    return result


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


def make_cluster_path(matrix_path, method):
    output_dir = os.path.dirname(matrix_path)
    file_base = os.path.splitext(os.path.basename(matrix_path))[0]
    

    if file_base.endswith("_matrix_lower"):
        file_base = file_base[:-13]
        
    return os.path.join(output_dir, f"{file_base}_cluster_{method}.csv")

def export(matrix_path, threshold, method):
    output_file = make_cluster_path(matrix_path, method)

    with open(matrix_path, "r") as temp_f:
        col_count = [len(l.split(",")) for l in temp_f.readlines()]
        column_names = [i for i in range(0, max(col_count))]

    df = pd.read_csv(
        matrix_path, delimiter=",", index_col=0, header=None, names=column_names
    )
    index = df.index.tolist()
    data = df.to_numpy()
    data = np.round(data, 2)

    df_result = get_clusters_dataframe(data, method, threshold, index)
    df_result.to_csv(output_file, index=False)

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
