import os
from collections import namedtuple
DocPaths = namedtuple(
    "DocPaths",
    [
        "stats",
        "full_matrix",
        "matrix_lower",
        "columns",
        "summary",
        "document",
        "cluster"
    ],
)

doc_dict = {
    "stats": "stats.csv",
    "full_matrix": "matrix.csv",
    "matrix_lower": "mat.csv",
    "columns": "columns.csv",
    "summary": "summary.csv",
    "document": "document.json",
    "cluster": "cluster.csv",   
}


def fetch_docpaths(tempdir_path):

    paths = {
        field: os.path.join(tempdir_path, filename)
        for field, filename in doc_dict.items()
    }
    return DocPaths(**paths)

