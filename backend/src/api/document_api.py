import os
import sys

sys.path.append(os.path.join(os.path.dirname(__file__), ".."))

from save_document import pack_document
from app_settings import add_recent_file
from document_state import save_document_settings
from utils import make_doc_id
from debug import open_doc_folder


def new_doc(new_document):
    id = make_doc_id()
    new_document(id)
    return id


def get_doc(doc_id: str, get_document):
    doc = get_document(doc_id)
    return doc._asdict() if doc else None


def save_doc(doc_id: str, path: str, save_as: bool, get_document, update_document):
    doc = get_document(doc_id)
    files = [
        os.path.join(doc.tempdir_path, file)
        for file in os.listdir(doc.tempdir_path)
    ]
    pack_document(path, files)
    if save_as == False:
        basename = os.path.basename(path)
        update_document(
            doc_id, filename=path, basename=basename, filetype="application/vnd.sdt"
        )
        add_recent_file(path)


def save_doc_settings(args: dict, get_document, update_document):
    doc = get_document(args["id"])
    update_document(
        args["id"],
        dataView=args["dataView"],
        heatmap=args["heatmap"],
        clustermap=args["clustermap"],
        distribution=args["distribution"],
    )
    doc = get_document(args["id"])
    save_document_settings(doc)


def close_doc(doc_id: str, remove_document):
    remove_document(doc_id)


def set_window_title(title: str, window, window_title_suffix):
    if window:
        window.title = f"{title}{window_title_suffix}"


def open_doc_folder_api(doc_id: str, get_document):
    doc = get_document(doc_id)
    return open_doc_folder(doc)
