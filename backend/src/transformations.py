from numpy import tril, triu_indices, around, nan
import numpy as np
import pandas as pd
from pandas import DataFrame


def to_triangle(matrix, convert_to_similarity=True, fill_value=nan):

    is_dataframe = isinstance(matrix, DataFrame)

    # Get numpy array and metadata
    if is_dataframe:
        data = matrix.to_numpy()
        index = matrix.index
        columns = matrix.columns
    else:
        data = matrix

    # Convert to similarity
    if convert_to_similarity:
        data = 100 - data

    # Create lower triangle
    result = tril(around(data, 2))

    # Fill upper triangle
    upper_indices = triu_indices(result.shape[0], k=1)
    if fill_value is None:
        result[upper_indices] = None
        result = np.where(np.isnan(result), None, result)
    else:
        result[upper_indices] = fill_value

    # Return in original format
    if is_dataframe:
        return DataFrame(result, index=index, columns=columns)
    else:
        return result


def similarity_triangle_to_matrix(matrix_lower: DataFrame) -> DataFrame:
    index = matrix_lower.index
    sim_lower = np.round(np.tril(100 - matrix_lower), 2)
    df = sim_lower + np.triu(sim_lower.T, 1)
    return DataFrame(df, index=index)


def read_csv_file(
    filepath, sep=",", header=None, index_col=None, as_list=False, extract_column=None
):
    # Special handling for matrix files with variable columns
    if header is None and index_col == 0 and sep == ",":
        with open(filepath, "r") as temp_f:
            col_count = [len(l.split(sep)) for l in temp_f.readlines()]
            column_names = [i for i in range(0, max(col_count))]
        df = pd.read_csv(
            filepath,
            delimiter=sep,
            index_col=index_col,
            header=header,
            names=column_names,
        )
    else:
        df = pd.read_csv(filepath, sep=sep, header=header, index_col=index_col)

    # Return based on options
    if extract_column:
        return df[extract_column].tolist()
    elif as_list:
        return df.values.tolist()
    else:
        return df


def read_csv_matrix(filepath):
    return read_csv_file(filepath, sep=",", header=None, index_col=0)


def read_stats_csv(filepath):
    return read_csv_file(filepath, header=0)


def read_columns_csv(filepath):
    # need typimg for id coliumn  **DtypeWarning: Columns (2) have mixed types. Specify dtype option on import or set low_memory=False.
    dtype_spec = {
        0: str,
        1: str,
        2: float,
    }  # First Sequence, Second Sequence, Identity Score
    df = pd.read_csv(filepath, header=0, dtype=dtype_spec)
    return df.values.tolist()


def read_tsv_file(filepath, id_column=None):
    return read_csv_file(filepath, sep="\t", extract_column=id_column)


def lzani_tsv_to_distance_matrix(results_tsv_path, ids_tsv_path, score_column="ani"):
    # Read IDs - explicitly get list
    ids_df = pd.read_csv(ids_tsv_path, sep="\t")
    all_ids = ids_df["id"].tolist()
    # Read results - explicitly get DataFrame
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
    
    # Set floor value to 0 for unreliable alignments
    # Values below 70% are considered unreliable according to the LZANI paper
    floor_value = 0
    
    # Replace unaligned NANs, 0s, AND values below 70% with 0
    matrix_np = np.where(np.isnan(matrix_np) | (matrix_np == 0) | (matrix_np < 20), floor_value, matrix_np)

    # LZANI has slightly different scores for query v. reference vs. reference v query. Make symmetric by averaging with transpose
    matrix_np = (matrix_np + matrix_np.T) / 2

    # fill diagonal with 100 for self-comparisons
    np.fill_diagonal(matrix_np, 100.0)

    # Convert similarity to distance
    distance_matrix = 100 - matrix_np

    return distance_matrix, all_ids
