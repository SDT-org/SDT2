from numpy import tril, triu_indices, around, nan
import numpy as np
from pandas.core.frame import DataFrame
from pandas import DataFrame

def dataframe_to_triangle(matrix: DataFrame) -> DataFrame:
    return matrix.where(np.tril(np.ones(matrix.shape, dtype=bool))).round(2)

def numpy_to_triangle(matrix: np.ndarray) -> np.ndarray:
    result = tril(around(matrix, 2))
    i_upper = np.triu_indices(result.shape[0], 1)
    result[i_upper] = None
    result = np.where(np.isnan(result), None, result)
    return result

def triangle_to_matrix(matrix_lower: DataFrame) -> DataFrame:
    # We need to round to 2 decimals for SDT 1 matrix files saving 100s as 100.001
    index = matrix_lower.index
    lower = np.round(np.tril(matrix_lower), 2)
    df = lower + np.triu(lower.T, 1)
    dataframe= DataFrame(df, index= index)
    return DataFrame(dataframe)
