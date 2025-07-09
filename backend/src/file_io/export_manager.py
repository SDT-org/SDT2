import shutil
import base64
import os
import json
from typing import List, Dict
from pandas.core.frame import DataFrame

from config.paths import (
    DATA_FILES,
    IMAGE_KEYS,
    DataKey,
    ImageFormat,
    ImageKey,
    build_document_paths,
)
from state.document_state import DocState

EXPORTABLE_DATA_KEYS: List[DataKey] = [
    "stats",
    "matrix",
    "triangle",
    "columns",
    "summary",
    "cluster_dir", 
]

MATRIX_ONLY_EXPORTABLE_DATA_KEYS: List[DataKey] = [
    "matrix",
    "columns",
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
        doc_path: str, export_path: str, prefix: str, image_format: ImageFormat, matrix_only = False
) -> List[List[str]]:
    paths_dict = build_document_paths(doc_path)._asdict()
    exportable_data_keys = (
        MATRIX_ONLY_EXPORTABLE_DATA_KEYS if matrix_only else EXPORTABLE_DATA_KEYS
    )
    data_file_paths = [
        [
            paths_dict[key],
            os.path.join(export_path, f"{prefix}_{DATA_FILES[key]}"),
        ]
        for key in exportable_data_keys
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
        if os.path.isdir(source):
            # If temp_target exists from a previous failed attempt, remove it before copytree
            if os.path.exists(temp_target):
                shutil.rmtree(temp_target)
            shutil.copytree(source, temp_target)
        else:
            shutil.copy2(source, temp_target)
        
        if os.path.isdir(target) and os.path.isdir(temp_target):
             shutil.rmtree(target) 
        elif os.path.isfile(target) and os.path.isdir(temp_target):
             os.remove(target) # Ensure target file is removed before replacing with a dir
        
        os.replace(temp_target, target)


def save_seq_dict_to_json(seq_dict, filename):
    with open(filename, "w") as file:
        json.dump(seq_dict, file, indent=4)


def save_run_settings_to_json(run_settings: Dict, path: str):
    with open(path, "w") as file:
        json.dump(run_settings, file, indent=4)
