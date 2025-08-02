import os
from app_settings import add_recent_file
from constants import default_window_title
from debug import open_doc_folder as debug_open_doc_folder
from document_state import save_document_settings
from save_document import pack_document
from utils import make_doc_id
from app_state import get_state, reset_state as app_reset_state, new_document, get_document, update_document, remove_document
from app_globals import assert_app_window


class Documents:
    def get_state_dict(self):
        return get_state()._asdict()

    def reset_state(self):
        app_reset_state()
        assert_app_window().title = default_window_title

    def new_doc(self):
        id = make_doc_id()
        new_document(id)
        return id

    def get_doc(self, doc_id: str):
        doc = get_document(doc_id)
        return doc._asdict() if doc else None

    def save_doc(self, doc_id: str, path: str, save_as: bool = False):
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

    def save_doc_settings(self, args: dict):
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

    def close_doc(self, doc_id: str):
        remove_document(doc_id)

    def set_window_title(self, title: str):
        from app import window_title_suffix
        assert_app_window().title = f"{title}{window_title_suffix}"

    def open_doc_folder(self, doc_id: str):
        doc = get_document(doc_id)
        return debug_open_doc_folder(doc)