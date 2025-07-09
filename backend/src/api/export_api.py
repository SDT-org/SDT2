import os
import sys
import urllib.parse
import mimetypes
import webview

sys.path.append(os.path.join(os.path.dirname(__file__), ".."))

from file_io.export_manager import (
    ImageFormat,
    ImageKey,
    build_source_target_pairs,
    do_export,
    save_image_from_api,
)
from document_paths import build_document_paths
from app_settings import update_app_settings
from utils import open_folder
from workflow import cluster


def save_svg_element(doc_id: str, selector: str, key: ImageKey, get_document, window):
    doc = get_document(doc_id)
    element = window.dom.get_element(selector)
    if element == None:
        raise Exception(f"could not find element: {selector}")
    inner_html = element.node["innerHTML"]
    data = f"<svg xmlns='http://www.w3.org/2000/svg'>{inner_html}</svg>"
    save_image_from_api(doc=doc, data=data, key=key, format="svg")
    return True


def save_svg_data(doc_id: str, data: str, key: ImageKey, format: ImageFormat, get_document):
    doc = get_document(doc_id)
    data = urllib.parse.unquote(data.split(",")[1])
    save_image_from_api(
        doc=doc,
        data=data,
        key=key,
        format=format,
    )
    return True


def save_raster_image(doc_id: str, data: str, key: ImageKey, format: ImageFormat, get_document):
    doc = get_document(doc_id)
    save_image_from_api(
        doc=doc,
        data=data,
        key=key,
        format=format,
    )
    return True


def export(args: dict, get_document, confirm_overwrite, matrix_filetypes):
    doc = get_document(args["doc_id"])
    update_app_settings(
        {
            "user_settings": {
                "export_path": args["export_path"],
                "open_folder_after_export": args["open_folder"],
            }
        }
    )
    doc_paths = build_document_paths(doc.tempdir_path)
    if args["output_cluster"] == True:
        cluster.export(
            matrix_path=doc_paths.matrix,
            cluster_data_output_dir=doc_paths.cluster_dir,
            seq_dict_path=doc_paths.seq_dict,
            threshold=args["cluster_threshold"],
            method=args["cluster_method"],
        )
    prefix_default = os.path.splitext(doc.basename)[0]
    args["prefix"] = args.get("prefix", prefix_default) or prefix_default
    filetype, _ = mimetypes.guess_type(doc.basename)
    matrix_only = filetype in matrix_filetypes
    source_target_pairs = build_source_target_pairs(
        doc.tempdir_path,
        args["export_path"],
        args["prefix"],
        args["image_format"],
        matrix_only,
    )
    source_paths, target_paths = zip(*source_target_pairs)
    for file in source_paths:
        if not os.path.exists(file):
            raise FileNotFoundError(f"File not found: {file}")
    if not confirm_overwrite(target_paths):
        return False
    do_export(source_target_pairs)
    if args["open_folder"]:
        open_folder(args["export_path"])
    return True


def confirm_overwrite(target_files, window):
    files_to_overwrite = [f for f in target_files if os.path.exists(f)]
    if not files_to_overwrite:
        return True
    message = (
        "The following files will be overwritten:\n\n"
        + "\n".join(map(os.path.basename, files_to_overwrite))
        + "\n\nAre you sure?"
    )
    return window.create_confirmation_dialog("Warning", message)
