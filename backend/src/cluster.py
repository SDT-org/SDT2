import os
import numpy as np
import pandas as pd
from collections import defaultdict
import networkx as nx

## switching from scipy to networkx. takes two threshold inputs now
def process_groups(data, index, threshold_1, threshold_2=0):
    # check for a threshold 2
    if threshold_2 is None or threshold_2 == 0:
        # set all values in the matrix that meet threshold 1 to 1 and all lower or NaN values to 0
        adjacency_matrix = (~np.isnan(data) & (data >= threshold_1)).astype(int)
        # Create a graph from the adjacency matrix
        G1 = nx.from_numpy_array(adjacency_matrix, parallel_edges=False, create_using=None)
        
        # empty dict to store the groups
        groups_dict = defaultdict(list)
        
        # look for connected components in graph G1
        for i, component in enumerate(nx.connected_components(G1)):
            for node_idx in component:
                groups_dict[i].append(index[node_idx])
        
        return groups_dict
    else:
        # Create two adjacency matrices, one for each threshold ~ is the bitwise 'not' operator
        adjacency_1 = (~np.isnan(data) & (data >= threshold_1)).astype(int)
        adjacency_2 = (~np.isnan(data) & (data >= threshold_2)).astype(int)
        
        # convert adjacency matrices to networkx graphs
        G1 = nx.from_numpy_array(adjacency_1)
        G2 = nx.from_numpy_array(adjacency_2)
        
        # find primary clusters with threshold_1
        groups_dict_1 = defaultdict(list)
        for i, component in enumerate(nx.connected_components(G1)):
            for node_idx in component:
                groups_dict_1[i].append(index[node_idx])
        
        # find subclusters with threshold_2 
        groups_dict_2 = defaultdict(list)
        for i, component in enumerate(nx.connected_components(G2)):
            for node_idx in component:
                groups_dict_2[i].append(index[node_idx])
        
        # return both cluster sets
        return groups_dict_1, groups_dict_2

def cluster_by_identity(clusters, nodes):
    output = []
    reverse_clusters = {}
    
    # create lookup table from sequence ID to primary cluster ID
    for group, values in clusters.items():
        for value in values:
            reverse_clusters[value] = group
    
    # initialize counters for subgroups within each primary cluster
    subgroup_counters = {group: 1 for group in clusters.keys()}
    
    # assign subgroups within each primary cluster
    for _, node_list in nodes.items():
        if node_list:
            # get first node to determine which primary cluster this belongs to
            first_value = node_list[0]
            if first_value in reverse_clusters:
                # get primary cluster ID (add 1 for human-readable indexing)
                group_number = reverse_clusters[first_value] + 1
                # get next available subgroup number for this primary cluster
                subgroup_number = subgroup_counters[reverse_clusters[first_value]]
                
                # process all nodes in this subcluster
                for value in node_list:
                    # only include if node belongs to the same primary cluster
                    if value in reverse_clusters:
                        output.append((value, group_number, subgroup_number))
                
                # increment subgroup counter for this primary cluster
                subgroup_counters[reverse_clusters[first_value]] += 1
    
    return output


def export(matrix_path, threshold_1=79, threshold_2=0, save_csv=True):
    output_dir = os.path.dirname(matrix_path)
    file_name = os.path.basename(matrix_path)
    file_base, _ = os.path.splitext(file_name)
    file_name = file_base.replace("_mat", "")
    output_file = os.path.join(output_dir, file_name + "_cluster.csv")
    
    with open(matrix_path, "r") as temp_f:
        col_count = [len(l.split(",")) for l in temp_f.readlines()]
        column_names = [i for i in range(0, max(col_count))]
    
    df = pd.read_csv(
        matrix_path, delimiter=",", index_col=0, header=None, names=column_names
    )
    
    index = df.index.tolist()
    data = df.to_numpy()
    data = np.round(data, 2)
    
    if threshold_2 != 0 and threshold_1 >= threshold_2:
        threshold_1, threshold_2 = threshold_2, threshold_1
    
    if threshold_2 is None or threshold_2 == 0:
        output = process_groups(data, index, threshold_1)
        flattened_output = [
            (item, key + 1) for key, sublist in output.items() for item in sublist
        ]
        df_result = pd.DataFrame(flattened_output)
        df_result.columns = ["SeqID", "Group - Threshold: " + str(threshold_1)]
    else:
        clusters, nodes = process_groups(data, index, threshold_1, threshold_2)
        output = cluster_by_identity(clusters, nodes)
        df_result = pd.DataFrame(output)
        df_result.columns = [
            "ID",
            "Group - Threshold: " + str(threshold_1),
            "Subgroup - Threshold: " + str(threshold_2),
        ]
    
    if save_csv:
        df_result.to_csv(output_file, index=False)
    
    return df_result

