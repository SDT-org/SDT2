from numpy import tril, triu_indices, around, nan
import numpy
from pandas.core.frame import DataFrame
from pandas import DataFrame

def dataframe_to_lower_triangle(matrix: DataFrame) -> DataFrame:
    dataframe = tril(around(matrix, 2))
    dataframe[triu_indices(dataframe.shape[0], k=1)] = nan

    return DataFrame(dataframe)

def numpy_to_lower_triangle(matrix: numpy.ndarray) -> numpy.ndarray:
    result = tril(around(matrix, 2))
    i_upper = numpy.triu_indices(result.shape[0], 1)
    result[i_upper] = None
    result = numpy.where(numpy.isnan(result), None, result)
    return result
