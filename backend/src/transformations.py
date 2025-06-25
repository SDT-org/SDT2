from numpy import tril, triu_indices, around, nan
import numpy as np
import pandas as pd
from pandas import DataFrame  


def dataframe_to_triangle(matrix: DataFrame) -> DataFrame:
    dataframe = tril(around(matrix, 2))
    dataframe[triu_indices(dataframe.shape[0], k=1)] = nan
    return DataFrame(dataframe)


def numpy_to_triangle(matrix: np.ndarray) -> np.ndarray:
    result = tril(around(matrix, 2))
    i_upper = np.triu_indices(result.shape[0], 1)
    result[i_upper] = None
    result = np.where(np.isnan(result), None, result)
    return result


def triangle_to_matrix(matrix_lower: DataFrame) -> DataFrame:
    """Convert a lower triangle matrix to a full symmetric matrix."""
    # We need to round to 2 decimals for SDT 1 matrix files saving 100s as 100.001
    index = matrix_lower.index
    lower = np.round(np.tril(matrix_lower), 2)
    df = lower + np.triu(lower.T, 1)
    dataframe = DataFrame(df, index=index)
    return dataframe

def lzani_tsv_to_distance_matrix(results_tsv_path, ids_tsv_path, score_column="ani"):
    """Convert LZANI TSV output files to a distance matrix."""
    # Read IDs and results
    all_ids = pd.read_csv(ids_tsv_path, sep="\t")["id"].tolist()
    df = pd.read_csv(results_tsv_path, sep="\t")
    
    # Create pivot table with similarity scores (still 0-1 scale)
    matrix = df.pivot_table(
        index="query", columns="reference", values=score_column, aggfunc="first"
    )
    
    # Reindex to ensure all IDs are present
    matrix = matrix.reindex(index=all_ids, columns=all_ids)
    
    # Convert to numpyand scale to percentage
    matrix_np = matrix.to_numpy() * 100
    
    # Replace unlaigned NANs or 0s with a biological floor to indicate no alignemnt
    matrix_np = np.where(np.isnan(matrix_np) | (matrix_np == 0), 25.0, matrix_np)
    
    # LZANI has sligntly difference scores for query v. reference vs. reverence v query..Make symmetric by averaging with transpose
    matrix_np = (matrix_np + matrix_np.T) / 2
    
    # fill diagonal with 100 for self-comparisons
    np.fill_diagonal(matrix_np, 100.0)
    
    # Convert similarity to distance
    distance_matrix = 100 - matrix_np
    
    return distance_matrix, all_ids
