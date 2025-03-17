import os
import numpy as np
from scipy.cluster.hierarchy import linkage, fcluster
import pandas as pd
from collections import defaultdict
from scipy.spatial.distance import squareform, pdist
from sklearn.manifold import MDS

def process_groups( data, index,threshold, method):
    # SET the upp triangle to zero
    i_upper = np.triu_indices(data.shape[0], 1)
    #transpost data to upper triangle
    data[i_upper] = data.T[i_upper]
    # convert to distance matrix
    distance_mat = 100 - data
    
    # create linkage matrix
    if method =="ward" or method == "centroid" or method == "complete":
        mds_coords = MDS(n_components=2, dissimilarity='precomputed', random_state=42).fit_transform(distance_mat)
        mds_distances = pdist(mds_coords)
        Z = linkage(mds_distances, method=method, metric='euclidean')
    else:
        condensed_dist = squareform(distance_mat)
        Z = linkage(condensed_dist, method=method, metric='precomputed')
        
    #set cut threshold
    cutby = 100 - threshold
    #identify clusters from threshold cut
    labels = fcluster(Z, t=cutby, criterion='distance')
    # store as dict
    groups_dict = defaultdict(list)
    for i, label in enumerate(labels):
        groups_dict[label-1].append(index[i])
    return groups_dict

def export(matrix_path, threshold, method, save_csv = True):
    output_dir = os.path.dirname(matrix_path)
    file_name = os.path.basename(matrix_path)
    file_base, _ = os.path.splitext(file_name)
    file_name = file_base.replace("_mat", "")
    output_file = os.path.join(output_dir, file_name + "_cluster_" + method + ".csv")
    with open(matrix_path, "r") as temp_f:
        col_count = [len(l.split(",")) for l in temp_f.readlines()]
        column_names = [i for i in range(0, max(col_count))]
    df = pd.read_csv(matrix_path, delimiter=",", index_col=0, header=None, names=column_names)
    index = df.index.tolist()
    data = df.to_numpy()
    data = np.round(data, 2) ## if we want to add a percision argument we cna do that here at some point
    output = process_groups(data, index,threshold, method)
    flattened_output = [(item, key + 1) for key, sublist in output.items() for item in sublist]
    df_result = pd.DataFrame(flattened_output)
    df_result.columns = ["ID", "Group - Threshold: " + str(threshold)]

    df_result.to_csv(output_file, index=False)
    return df_result
