import json
import csv
import os
from typing import Dict
import numpy as np
from pandas.core.frame import DataFrame
from transformations import to_triangle
from workflow.models import RunSettings


def save_matrix_to_csv(df, matrix_path, triangle_path):
    print(f"Exporting matrix files for {len(df)} sequences...")
    # to_triangle returns a DataFrame when given a DataFrame
    triangle = to_triangle(df)
    # stupid type error
    if isinstance(triangle, DataFrame):
        print(f"Writing triangle matrix: {os.path.basename(triangle_path)}")
        triangle.to_csv(triangle_path, mode="wt", header=False, index=True, sep=",")
    print(f"Writing full matrix: {os.path.basename(matrix_path)}")
    df.to_csv(matrix_path, mode="w", header=False, index=True, sep=",")
    print("Matrix export completed")


# takes dissimilarity matrix and converts it to a columnar format
def save_cols_to_csv(df, path):
    print(f"Exporting columnar data for {len(df)} sequences to {os.path.basename(path)}...")
    order = df.index
    df.columns = df.index
    
    # Convert DataFrame to numpy array for faster access
    matrix = df.to_numpy()
    n = len(order)
    
    # Calculate expected number of pairs
    num_pairs = n * (n - 1) // 2
    print(f"Generating {num_pairs} pairwise comparisons...")
    
    # Get indices for lower triangle (i > j)
    i_indices, j_indices = np.tril_indices(n, k=-1)
    
    # Extract values using vectorized operations
    print("Extracting similarity scores...")
    row_ids = [order[i] for i in i_indices]
    col_ids = [order[j] for j in j_indices]
    # Convert distance to similarity: 100 - distance
    identity_scores = 100 - matrix[i_indices, j_indices]
    
    # Write directly to CSV for better performance
    print(f"Writing CSV file with {len(row_ids)} rows...")
    with open(path, 'w', newline='', encoding='utf-8') as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(["First Sequence", "Second Sequence", "Identity Score"])
        writer.writerows(zip(row_ids, col_ids, identity_scores))
    
    print(f"Columnar export completed: {os.path.basename(path)}")


def save_stats_to_csv(seq_stats, filename):
    print(f"Exporting sequence statistics for {len(seq_stats)} sequences to {os.path.basename(filename)}...")
    stats_list = []
    for key, value in seq_stats.items():
        stats_list.append([key, value[0], value[1]])
    stats_df = DataFrame(stats_list, columns=["Sequence", "GC %", "Sequence Length"])
    stats_df.to_csv(filename, mode="w", header=True, index=False, sep=",")
    print("Statistics export completed")


def save_seq_dict_to_json(seq_dict, filename):
    print(f"Exporting sequence dictionary for {len(seq_dict)} sequences to {os.path.basename(filename)}...")
    with open(filename, "w") as file:
        json.dump(seq_dict, file, indent=4)
    print("Sequence dictionary export completed")


def save_run_settings_to_json(run_settings: Dict, path: str):
    with open(path, "w") as file:
        json.dump(run_settings, file, indent=4)
