from pandas.core.frame import DataFrame
from .transformations import to_triangle


def save_matrix_to_csv(df, matrix_path, triangle_path):
    # to_triangle returns a DataFrame when given a DataFrame
    triangle = to_triangle(df)
    # stupid type error
    if isinstance(triangle, DataFrame):
        triangle.to_csv(triangle_path, mode="wt", header=False, index=True, sep=",")
    df.to_csv(matrix_path, mode="w", header=False, index=True, sep=",")


# takes dissimilarity matrix and converts it to a columnar format
def save_cols_to_csv(df, path):
    order = df.index
    df.columns = df.index
    columnar_output = []
    for i_idx, row_id in enumerate(order):
        for j_idx, col_id in enumerate(order):
            if i_idx > j_idx:
                # convert back to a sim matrix - distance to similarity
                columnar_output.append([row_id, col_id, 100 - df.iloc[i_idx, j_idx]])
    columnar_df = DataFrame(
        columnar_output, columns=["First Sequence", "Second Sequence", "Identity Score"]
    )
    columnar_df.to_csv(path, mode="w", header=True, index=False, sep=",")


def save_stats_to_csv(seq_stats, filename):
    stats_list = []
    for key, value in seq_stats.items():
        stats_list.append([key, value[0], value[1]])
    stats_df = DataFrame(stats_list, columns=["Sequence", "GC %", "Sequence Length"])
    stats_df.to_csv(filename, mode="w", header=True, index=False, sep=",")
