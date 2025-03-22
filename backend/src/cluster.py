import os
import numpy as np
from scipy.cluster.hierarchy import linkage, fcluster,dendrogram
import pandas as pd
from collections import defaultdict
from scipy.spatial.distance import squareform, pdist
from sklearn.manifold import MDS

## creata blank cache which is dict of dicts to store hash ids with a dict of matrices
linkage_cache = {}
def calculate_linkages(data, index):

    global linkage_cache

    methods = [ "single", "complete", "average", "weighted", "centroid", "median","ward"]
    ## create uniquer id for each linkage cache
    id_hash= (hash(tuple(index)))

    # check if id in linkage dict, if not add key with blank val
    if id_hash not in linkage_cache:
        linkage_cache[id_hash] = {}

    # default check to false...assume linkage avail
    run_linkage=False


    #loop throuhgh each method if not in cachec add id to dict with blank key, or break loop
    for method in methods:
        if method not in linkage_cache[id_hash]:
            run_linkage=True
            break
    # if no linkage is needed, empty retrun as we want linkage_cache wiand it should be populated now
    if not run_linkage:
        return

    # SET the upp triangle to zero
    i_upper = np.triu_indices(data.shape[0], 1)
    #transpost data to upper triangle
    data[i_upper] = data.T[i_upper]
    # convert to distance matrix
    distance_mat = 100 - data

     #Redduce dimensionality with MDS once and store for neccessary linkages
    mds_coords = MDS(n_components=2, dissimilarity='precomputed', random_state=42).fit_transform(distance_mat)  # 42 is the answer to everything
    mds_distances = pdist(mds_coords)

    #convert to condensed distance matrix for straight likange
    condensed_dist = squareform(distance_mat)

    #Pass through each method and calculate linkage if not already in cache
    for method in methods:
        #check if id is in cache
        if method in linkage_cache[id_hash]:
             continue
        # create linkage matrix
        if method =="ward" or method == "centroid" or method == "complete":
            Z = linkage(mds_distances, method=method, metric='euclidean')
        else:
            Z = linkage(condensed_dist, method=method, metric='precomputed')

        linkage_cache[id_hash][method] = Z


def get_linkage(index, method):
    #recreate hash id to match cache
    id_hash = hash(tuple(index))
    #check if id in cache
    if id_hash in linkage_cache and method in linkage_cache[id_hash]:
        print(f"Using cached {method} for dataset {id_hash}")
        return linkage_cache[id_hash][method]
    else:
        print(f"ID not found in cacche for {method} linkage for dataset {id_hash}")
        return None

def get_linkage_method_order(data, method, index):
    Z = get_linkage_matrix(data, method)
    dendro = dendrogram(Z, no_plot=True)
    leaf_indices = dendro['leaves']
    new_order = [index[i] for i in leaf_indices]

    return new_order

def get_linkage_matrix(data, method):
    print(data, method)
    # # SET the upp triangle to zero
    # i_upper = np.triu_indices(data.shape[0], 1)
    # #transpost data to upper triangle
    # data[i_upper] = data.T[i_upper]
    # # convert to distance matrix
    # distance_mat = 100 - data

     #Redduce dimensionality with MDS once and store for neccessary linkages
    mds_coords = MDS(n_components=2, dissimilarity='precomputed', random_state=42).fit_transform(data)  # 42 is the answer to everything
    mds_distances = pdist(mds_coords)

    #convert to condensed distance matrix for straight likange
    condensed_dist = squareform(data)

    if method =="ward" or method == "centroid" or method == "complete":
        Z = linkage(mds_distances, method=method, metric='euclidean')
    else:
        Z = linkage(condensed_dist, method=method, metric='precomputed')

    return Z

def get_cluster_data(Z, threshold, index):
    #set cut threshold
    cutby = 100 - threshold
    #identify clusters from threshold cut
    clusters = fcluster(Z, t=cutby, criterion='distance')
    # store as dict
    cluster_dict = defaultdict(list)
    for i, label in enumerate(clusters):
        cluster_dict[label-1].append(index[i])

    return cluster_dict

def export(matrix_path, threshold, method, save_csv = True):

   #set variables
    output_dir = os.path.dirname(matrix_path)
    file_name = os.path.basename(matrix_path)
    file_base, _ = os.path.splitext(file_name)
    file_name = file_base.replace("_mat", "")


    #parse data
    output_file = os.path.join(output_dir, file_name + "_cluster_" + method + ".csv")
    with open(matrix_path, "r") as temp_f:
        col_count = [len(l.split(",")) for l in temp_f.readlines()]
        column_names = [i for i in range(0, max(col_count))]

    df = pd.read_csv(matrix_path, delimiter=",", index_col=0, header=None, names=column_names)
    index = df.index.tolist()
    data = df.to_numpy()
    data = np.round(data, 2) ## if we want to add a percision argument we cna do that here at some point
    Z = get_linkage(index, method)

    # If not cached, calculate all linkages
    if Z is None:
        calculate_linkages(data, index)
        Z = get_linkage(index, method)

    output = get_cluster_data(Z, threshold, index)
    flattened_output = [(item, key + 1) for key, sublist in output.items() for item in sublist]
    df_result = pd.DataFrame(flattened_output)
    df_result.columns = ["ID", "Group - Threshold: " + str(threshold)]

    df_result.to_csv(output_file, index=False)
    return df_result
