from numpy import tril, triu_indices, around, nan
import numpy as np
import pandas as pd
from pandas import DataFrame


def dataframe_to_triangle(matrix: DataFrame) -> DataFrame:
    # Convert distance to similarity
    similarity_matrix = 100 - matrix
    
    # Create lower triangle
    dataframe = tril(around(similarity_matrix, 2))
    dataframe[triu_indices(dataframe.shape[0], k=1)] = nan
    
    return DataFrame(dataframe, index=matrix.index, columns=matrix.columns)


def numpy_to_triangle(matrix: np.ndarray) -> np.ndarray:
    # Convert distance to similarity
    similarity_matrix = 100 - matrix
    
    # Create lower triangle
    result = tril(around(similarity_matrix, 2))
    i_upper = np.triu_indices(result.shape[0], 1)
    result[i_upper] = None
    result = np.where(np.isnan(result), None, result)
    return result


def triangle_to_matrix(matrix_lower: DataFrame) -> DataFrame:
    index = matrix_lower.index
    lower = np.round(np.tril(matrix_lower), 2)
    df = lower + np.triu(lower.T, 1)
    dataframe = DataFrame(df, index=index)
    return dataframe

def read_csv_matrix(filepath):
    with open(filepath, "r") as temp_f:
        col_count = [len(l.split(",")) for l in temp_f.readlines()]
        column_names = [i for i in range(0, max(col_count))]
    
    return pd.read_csv(
        filepath, delimiter=",", index_col=0, header=None, names=column_names
    )


def lzani_tsv_to_distance_matrix(results_tsv_path, ids_tsv_path, score_column="ani"):
    all_ids = pd.read_csv(ids_tsv_path, sep="\t")["id"].tolist()
    df = pd.read_csv(results_tsv_path, sep="\t")
    
    # Create pivot table with similarity scores (still 0-1 scale)
    if df.empty:
        matrix = pd.DataFrame(index=all_ids, columns=all_ids)
    else:
        matrix = df.pivot_table(
            index="query", columns="reference", values=score_column, aggfunc="first"
        )
    
    # Reindex to ensure all IDs are present
    matrix = matrix.reindex(index=all_ids, columns=all_ids)
    
    # Convert to numpy and scale to percentage
    matrix_np = matrix.to_numpy() * 100
    
    # Find the minimum non-zero, non-NaN value to use as floor
    valid_values = matrix_np[~np.isnan(matrix_np) & (matrix_np > 0)]
    if len(valid_values) > 0:
        min_value = np.min(valid_values)
        # Use 90% of minimum value as floor for unaligned sequences
        floor_value = min_value * 0.9
    else:
        # Fallback if no valid values found
        floor_value = 25.0
    
    # Replace unaligned NANs or 0s with the calculated floor
    matrix_np = np.where(np.isnan(matrix_np) | (matrix_np == 0), floor_value, matrix_np)
    
    # LZANI has slightly different scores for query v. reference vs. reference v query. Make symmetric by averaging with transpose
    matrix_np = (matrix_np + matrix_np.T) / 2
    
    # fill diagonal with 100 for self-comparisons
    np.fill_diagonal(matrix_np, 100.0)
    
    # Convert similarity to distance
    distance_matrix = 100 - matrix_np
    
    return distance_matrix, all_ids
