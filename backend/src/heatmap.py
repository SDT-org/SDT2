from numpy import tril, triu_indices, around, nan
import numpy as np
from pandas.core.frame import DataFrame
from pandas import DataFrame

def dataframe_to_lower_triangle(matrix: DataFrame) -> DataFrame:
    dataframe = tril(around(matrix, 2))
    dataframe[triu_indices(dataframe.shape[0], k=1)] = nan

    return DataFrame(dataframe)

def numpy_to_lower_triangle(matrix: np.ndarray) -> np.ndarray:
    result = tril(around(matrix, 2))
    i_upper = np.triu_indices(result.shape[0], 1)
    result[i_upper] = None
    result = np.where(np.isnan(result), None, result)
    return result

def lower_triangle_to_full_matrix(lower_triangle: DataFrame) -> DataFrame:
    lower = np.tril(lower_triangle)
    dataframe = lower + np.triu(lower.T, 1)

    return DataFrame(dataframe)
