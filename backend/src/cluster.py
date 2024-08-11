import os
import numpy as np
from scipy.sparse import csr_matrix
from scipy.sparse.csgraph import connected_components
import pandas as pd
from collections import defaultdict

def process_groups(threshold, data, index):
    # create adjacensy matrix to id which cells are related by the threshold marking as binary with (1) for related or (0) for not meeting the threshold
    adjacency_matrix = (data >= threshold).astype(int)
    # create sparse matrix(absent 0s) for memeory efficiancy
    sparse_matrix = csr_matrix(adjacency_matrix)
    # identify connected components
    _, labels = connected_components(
        csgraph=sparse_matrix, directed=False, return_labels=True
    )
    groups_dict = defaultdict(list)
    for i, label in enumerate(labels):
        groups_dict[label].append(index[i])
    groups = {}
    for indx, clade in enumerate(labels):
        groups.update({indx: clade})
    return groups_dict


def cluster_by_identity(clusters, nodes):
    output = []
    reverse_clusters = {}

    # Create reverse lookup dictionary were values in value list are extracted to key and groups are assigned to value
    for group, values in clusters.items():
        for value in values:
            reverse_clusters[value] = group

    # Initialize subgroup counters
    subgroup_counters = {group: 1 for group in clusters.keys()}

    # Iterate through nodes to determine the subgroup_number within each group_number
    for node_key, node_list in nodes.items():
        if node_list:
            first_value = node_list[0]
            if first_value in reverse_clusters:
                group_number = reverse_clusters[first_value] + 1
                subgroup_number = subgroup_counters[reverse_clusters[first_value]]
                for value in node_list:
                    if value in reverse_clusters:
                        output.append((value, group_number, subgroup_number))
                subgroup_counters[reverse_clusters[first_value]] += 1

    return output



def export(matrix_path, threshold_1=79, threshold_2=0):
    output_dir = os.path.dirname(matrix_path)
    file_name = os.path.basename(matrix_path)
    file_base, _ = os.path.splitext(file_name)
    file_name = file_base.replace("_mat", "")
    output_file = os.path.join(output_dir, file_name + "_cluster.csv")

    # https://stackoverflow.com/a/57824142
    # SDT1 matrix CSVs do not have padding for columns
    with open(matrix_path, 'r') as temp_f:
        col_count = [ len(l.split(",")) for l in temp_f.readlines() ]
        column_names = [i for i in range(0, max(col_count))]

    df = pd.read_csv(matrix_path, delimiter=",", index_col=0, header=None, names=column_names)
    # extract index
    index = df.index.tolist()
    # convert df to np array
    data = df.to_numpy()
    # format values
    data = np.round(data, 2)
    # maintain order of threshold processing
    if threshold_2 != 0 and threshold_1 >= threshold_2:
        threshold_1, threshold_2 = threshold_2, threshold_1
    # handle instances of no threshold_2
    if threshold_2 is None or threshold_2 == 0:
        output = process_groups(threshold_1, data, index)
        flattened_output = [
            (item, key + 1) for key, sublist in output.items() for item in sublist
        ]
        df = pd.DataFrame(flattened_output)
        df.columns = ["ID", "Group 1 - Theshold: " + str(threshold_1)]

    else:
        clusters = process_groups(threshold_1, data, index)
        nodes = process_groups(threshold_2, data, index)
        output = cluster_by_identity(clusters, nodes)
        df = pd.DataFrame(output)
        df.columns = [
            "ID",
            "Group 1 - Theshold: " + str(threshold_1),
            "Group 2 - Theshold: " + str(threshold_2),
        ]
    df.to_csv(output_file, index=False)
