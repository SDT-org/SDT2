import os
from collections import namedtuple
from typing import Literal, Dict, List

DataKey = Literal[
    "stats", "matrix", "triangle", "columns", "summary", "settings", "cluster"
]

ImageKey = Literal["heatmap", "clustermap", "distribution_histogram", "distribution_violin"]
ImageFormat = Literal["png", "jpeg", "svg"]

DATA_FILES: Dict[DataKey, str] = {
    "stats": "stats.csv",
    "matrix": "matrix.csv",
    "triangle": "triangle.csv",
    "columns": "columns.csv",
    "summary": "summary.txt",
    "settings": "settings.json",
    "cluster": "cluster.csv",
}

IMAGE_KEYS: List[ImageKey] = ["heatmap", "clustermap", "distribution_histogram", "distribution_violin"]
IMAGE_FORMATS: List[ImageFormat] = ["png", "jpeg", "svg"]

ImageFormatPaths = namedtuple("ImageFormatPaths", IMAGE_FORMATS)
ImagePaths = namedtuple("ImagePaths", IMAGE_KEYS)
DocumentPaths = namedtuple("DocumentPaths", list(DATA_FILES.keys()) + ["images"])


def build_document_paths(doc_path: str) -> DocumentPaths:
    data_file_paths = {
        field: os.path.join(doc_path, filename)
        for field, filename in DATA_FILES.items()
    }

    image_paths = {}
    for key in IMAGE_KEYS:
        format_paths = {}
        for format in IMAGE_FORMATS:
            format_paths[format] = os.path.join(doc_path, f"{key}.{format}")
        image_paths[key] = ImageFormatPaths(**format_paths)

    return DocumentPaths(**data_file_paths, images=ImagePaths(**image_paths))
