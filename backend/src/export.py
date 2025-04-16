import shutil
from typing import List
import base64
import os
from document_paths import (
    DATA_FILES,
    IMAGE_KEYS,
    DataKey,
    ImageFormat,
    ImageKey,
    build_document_paths,
)
from document_state import DocState

EXPORTABLE_DATA_KEYS: List[DataKey] = [
    "stats",
    "matrix",
    "triangle",
    "columns",
    "summary",
    "cluster",
]


def save_image_from_api(doc: DocState, data: str, key: ImageKey, format: ImageFormat):
    target = (
        build_document_paths(doc.tempdir_path).images._asdict()[key]._asdict()[format]
    )

    if format == "svg":
        with open(target, "w", encoding="utf-8") as file:
            file.write(data)
    else:
        bytes = base64.b64decode(data.split(",")[1])
        with open(target, "wb") as file:
            file.write(bytes)


def build_source_target_pairs(
    doc_path: str, export_path: str, prefix: str, image_format: ImageFormat
) -> List[List[str]]:
    paths_dict = build_document_paths(doc_path)._asdict()
    data_file_paths = [
        [
            paths_dict[key],
            os.path.join(export_path, f"{prefix}_{DATA_FILES[key]}"),
        ]
        for key in EXPORTABLE_DATA_KEYS
    ]
    image_paths = [
        [
            getattr(paths_dict["images"], key)._asdict()[image_format],
            os.path.join(export_path, f"{prefix}_{key}.{image_format}"),
        ]
        for key in IMAGE_KEYS
    ]

    return data_file_paths + image_paths


def do_export(source_target_pairs):
    for source, target in source_target_pairs:
        temp_target = target + ".tmp"
        shutil.copy2(source, temp_target)
        os.replace(temp_target, target)
