import os
import numpy as np
from scipy.cluster.hierarchy import linkage, fcluster, dendrogram
import pandas as pd
from collections import defaultdict
from scipy.spatial.distance import squareform, pdist
from sklearn.manifold import MDS
import time

# Global cache for linkage matrices
# Structure: {hash_id: {'methods': {method: Z_matrix}, 'timestamp': last_accessed_time}}
linkage_cache = {}

def calculate_linkages(data, index):
    global linkage_cache
    
    methods = ["single", "complete", "average", "weighted", "centroid", "median", "ward"]
    id_hash = hash(tuple(index))
    
    # Update access timestamp or initialize cache entry
    current_time = time.time()
    
    if id_hash not in linkage_cache:
        linkage_cache[id_hash] = {'methods': {}, 'timestamp': current_time}
    else:
        linkage_cache[id_hash]['timestamp'] = current_time
    
    # Check if all methods are already in cache
    all_cached = True
    for method in methods:
        if method not in linkage_cache[id_hash]['methods']:
            all_cached = False
            break
    
    # Early return if all methods are cached
    if all_cached:
        return
    
    # Prepare the distance matrix
    dist_data = data.copy()
    
    # Set the upper triangle to mirror lower triangle
    i_upper = np.triu_indices(dist_data.shape[0], 1)
    dist_data[i_upper] = dist_data.T[i_upper]
    
    # Check if it's a similarity matrix (values near 100) or distance matrix
    if np.max(dist_data) > 90:
        distance_mat = 100 - dist_data
    else:
        distance_mat = dist_data
    
    # MDS for methods that need Euclidean distance
    mds_coords = MDS(n_components=2, dissimilarity='precomputed', random_state=42).fit_transform(distance_mat)
    mds_distances = pdist(mds_coords)
    
    # Convert to condensed distance matrix
    condensed_dist = squareform(distance_mat)
    
    # Calculate linkage for each method not in cache
    for method in methods:
        if method in linkage_cache[id_hash]['methods']:
            continue
            
        if method in ["ward", "centroid", "median"]:
            Z = linkage(mds_distances, method=method, metric='euclidean')
        else:
            Z = linkage(condensed_dist, method=method, metric='precomputed')
        
        linkage_cache[id_hash]['methods'][method] = Z

def get_linkage(index, method):
    global linkage_cache
    
    id_hash = hash(tuple(index))
    
    if id_hash in linkage_cache and method in linkage_cache[id_hash]['methods']:
        # Update timestamp when accessed
        linkage_cache[id_hash]['timestamp'] = time.time()
        return linkage_cache[id_hash]['methods'][method]
    else:
        return None

def get_linkage_method_order(data, method, index):
    # Make sure all linkage methods are calculated and cached
    calculate_linkages(data, index)
    
    # Get hash ID for the cache lookup
    id_hash = hash(tuple(index))
    
    # Get linkage matrix for the specified method
    if id_hash in linkage_cache and method in linkage_cache[id_hash]['methods']:
        Z = linkage_cache[id_hash]['methods'][method]
    else:
        print(f"Warning: Linkage matrix for {method} not found in cache")
        return index  # Return original order as fallback
    
    # Create dendrogram to get leaf order
    dendro = dendrogram(Z, no_plot=True)
    leaf_indices = dendro['leaves']
    new_order = [index[i] for i in leaf_indices]
    
    return new_order

def purge_old_cache_entries(max_age_seconds=3600, max_entries=10):
    global linkage_cache
    
    current_time = time.time()
    
    # Remove entries older than max_age_seconds
    old_keys = [k for k, v in linkage_cache.items() 
                if current_time - v['timestamp'] > max_age_seconds]
    
    for k in old_keys:
        del linkage_cache[k]
    
    # If still too many entries, remove oldest ones
    if len(linkage_cache) > max_entries:
        sorted_keys = sorted(linkage_cache.keys(), 
                             key=lambda k: linkage_cache[k]['timestamp'])
        for k in sorted_keys[:len(linkage_cache) - max_entries]:
            del linkage_cache[k]

def export(matrix_path, threshold, method, save_csv=True):
    # Set variables
    output_dir = os.path.dirname(matrix_path)
    file_name = os.path.basename(matrix_path)
    file_base, _ = os.path.splitext(file_name)
    file_name = file_base.replace("_mat", "")

    # Parse data
    output_file = os.path.join(output_dir, file_name + "_cluster_" + method + ".csv")
    
    with open(matrix_path, "r") as temp_f:
        col_count = [len(l.split(",")) for l in temp_f.readlines()]
        column_names = [i for i in range(0, max(col_count))]

    df = pd.read_csv(matrix_path, delimiter=",", index_col=0, header=None, names=column_names)
    index = df.index.tolist()
    data = df.to_numpy()
    data = np.round(data, 2)
    
    # Check if the linkage is already in the cache
    id_hash = hash(tuple(index))
    
    # If not in cache, calculate all linkages
    if id_hash not in linkage_cache or method not in linkage_cache.get(id_hash, {}).get('methods', {}):
        calculate_linkages(data, index)
    
    # Get the linkage matrix
    Z = linkage_cache[id_hash]['methods'][method]
    
    # Get cluster data
    output = get_cluster_data(Z, threshold, index)
    flattened_output = [(item, key + 1) for key, sublist in output.items() for item in sublist]
    df_result = pd.DataFrame(flattened_output)
    df_result.columns = ["ID", "Group - Threshold: " + str(threshold)]

    if save_csv:
        df_result.to_csv(output_file, index=False)
    
    # Periodically clean up old cache entries
    purge_old_cache_entries()
    
    return df_result

def get_cluster_data(Z, threshold, index):
    # Set cut threshold
    cutby = 100 - threshold
    # Identify clusters from threshold cut
    clusters = fcluster(Z, t=cutby, criterion='distance')
    # Store as dict
    cluster_dict = defaultdict(list)
    for i, label in enumerate(clusters):
        cluster_dict[label-1].append(index[i])
    
    return cluster_dict

def reorder_clusters(values, clusters):
    # Let's make it NP
    labels = np.array(clusters)
    # Create newly reordered labels that make sense
    new_order_labels = np.searchsorted(np.unique(labels), labels)
    return (new_order_labels + 1).tolist()