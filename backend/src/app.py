import os
import sys
import platform
import mimetypes
import json
import urllib.parse
from tempfile import TemporaryDirectory
from typing import Dict
import webview
import psutil
import numpy as np
from pandas import DataFrame
from numpy import eye, where, nan, nanmin, nanmax
from multiprocessing import Manager, cpu_count as get_cpu_count
from shutil import copy

current_file_path = os.path.abspath(os.path.join(os.path.dirname(__file__)))
sys.path.append(os.path.join(current_file_path, "../../"))
sys.path.append(os.path.join(current_file_path, "."))

from workflow.models import WorkflowResult, WorkflowRun, RunSettings, LzaniSettings, ParasailSettings
from app_state import create_app_state
from constants import default_window_title, matrix_filetypes
from utils import make_doc_id, get_child_process_info, open_folder
from api.webview_windows import get_html_path, about_window, manual_window
from app_settings import (
    add_recent_file,
    load_app_settings,
    remove_recent_file,
    save_app_settings,
    update_app_settings,
)
from save_document import pack_document, unpack_document
from document_state import save_document_settings
from document_paths import ImageKey, build_document_paths
from transformations import (
    to_triangle,
    similarity_triangle_to_matrix,
    read_csv_matrix,
    read_stats_csv,
)
from export_utils import save_cols_to_csv
from export import (
    ImageFormat,
    build_source_target_pairs,
    do_export,
    save_image_from_api,
)
from file_utils import read_json_file
from debug import open_doc_folder
from workflow.runner import run_process, run_parse
from workflow import cluster
from platformdirs import user_documents_dir
from config import app_version, dev_frontend_host

# Import API modules
from api import file_operations, workflow_api, document_api, export_api, data_api, system_api

is_compiled = "__compiled__" in globals()
is_macos = platform.system() == "Darwin"
is_windows = platform.system() == "Windows"
temp_dir = TemporaryDirectory()
window_title_suffix = "" if is_macos else " - SDT2"
try:
    cpu_count = get_cpu_count()
except:
    cpu_count = 1
mimetypes.add_type("text/fasta", ".fasta")
mimetypes.add_type("text/fasta", ".fas")
mimetypes.add_type("text/fasta", ".faa")
mimetypes.add_type("text/fasta", ".fnt")
mimetypes.add_type("text/fasta", ".fa")
mimetypes.add_type("application/vnd.sdt", ".sdt")
window = None
canceled = None
start_time = None  # TODO: move into workflow
workflow_runs: Dict[str, WorkflowRun] = {}
parsed_workflow_results: Dict[str, WorkflowResult] = {}


def assert_window():
    assert window is not None
    return window




class Api:
    def close_app(self):
        return system_api.close_app()

    def show_about(self):
        return system_api.show_about()

    def show_manual(self):
        return system_api.show_manual()

    def app_config(self):
        return system_api.app_config()

    def app_settings(self):
        return system_api.app_settings()

    def processes_info(self):
        return system_api.processes_info()

    def get_available_memory(self):
        return system_api.get_available_memory()

    def open_file(self, filepath: str, doc_id: str | None = None):
        return file_operations.handle_open_file(
            filepath, doc_id, temp_dir, get_state, new_document, 
            make_doc_id, parsed_workflow_results, remove_empty_documents
        )

    def open_file_dialog(self, doc_id: str | None = None, filepath: str | None = None):
        return file_operations.open_file_dialog(
            doc_id, filepath, lambda fp, di: self.open_file(fp, di)
        )

    def save_file_dialog(self, filename: str, directory: str | None = None):
        return file_operations.save_file_dialog(filename, directory)

    def select_path_dialog(self, directory: str | None = None):
        return file_operations.select_path_dialog(directory)

    def get_state(self):
        return system_api.get_state(get_state)

    def reset_state(self):
        return system_api.reset_state(reset_state, assert_window(), default_window_title)

    def start_workflow_run(self, args: dict):
        return workflow_api.start_workflow_run(
            args, get_state, get_document, update_document, 
            set_state, parsed_workflow_results, workflow_runs
        )

    def get_workflow_run_status(self, doc_id: str):
        return workflow_api.get_workflow_run_status(doc_id, workflow_runs)

    def cancel_run(self, doc_id: str, run_settings: str):
        global canceled
        return workflow_api.cancel_run(
            doc_id, run_settings, canceled, get_state, 
            update_document, set_state, workflow_runs, reset_state, assert_window()
        )

    def confirm_overwrite(self, target_files):
        return export_api.confirm_overwrite(target_files, assert_window())

    def save_svg_element(self, doc_id: str, selector: str, key: ImageKey):
        return export_api.save_svg_element(doc_id, selector, key, get_document, assert_window())

    def save_svg_data(self, doc_id: str, data: str, key: ImageKey, format: ImageFormat):
        return export_api.save_svg_data(doc_id, data, key, format, get_document)

    def save_raster_image(
        self, doc_id: str, data: str, key: ImageKey, format: ImageFormat
    ):
        return export_api.save_raster_image(doc_id, data, key, format, get_document)

    def export(self, args: dict):
        return export_api.export(
            args, get_document, lambda files: self.confirm_overwrite(files), matrix_filetypes
        )

    def get_data(self, doc_id: str):
        return data_api.get_data(doc_id, get_document, update_document)

    def get_clustermap_data(self, doc_id: str, threshold: float, method: str):
        return data_api.get_clustermap_data(doc_id, threshold, method, get_document)

    def new_doc(self):
        return document_api.new_doc(new_document)

    def get_doc(self, doc_id: str):
        return document_api.get_doc(doc_id, get_document)

    def save_doc(self, doc_id: str, path: str, save_as: bool = False):
        return document_api.save_doc(doc_id, path, save_as, get_document, update_document)

    def save_doc_settings(self, args: dict):
        return document_api.save_doc_settings(args, get_document, update_document)

    def close_doc(self, doc_id: str):
        return document_api.close_doc(doc_id, remove_document)

    def set_window_title(self, title: str):
        return document_api.set_window_title(title, assert_window(), window_title_suffix)

    def open_doc_folder(self, doc_id: str):
        return document_api.open_doc_folder_api(doc_id, get_document)


def push_backend_state(window: webview.Window):
    if window is None:
        # This can happen if the state updates before the window is fully initialized
        print(
            "Warning: push_backend_state called with window=None. Skipping UI update."
        )
        return
    state = get_state()
    dict_state = lambda t: {
        f: [d._asdict() for d in getattr(t, f)] if f == "documents" else getattr(t, f)
        for f in t._fields
    }
    js_app_state = json.dumps(dict(state=dict_state(state)))
    window.evaluate_js(
        f"document.dispatchEvent(new CustomEvent('sync-state', {{ detail: {js_app_state} }}))"
    )


def on_closed():
    # map(lambda doc: doc.cleanup(), get_state().documents)
    global canceled
    workflow_api.do_cancel_run(canceled, get_state, update_document, set_state, workflow_runs)
    os._exit(0)


if __name__ == "__main__":
    api = Api()
    window = webview.create_window(
        default_window_title,
        url=get_html_path(),
        js_api=api,
        # TODO: store last window size and position
        width=1200,
        height=900,
        min_size=(900, 800),
        confirm_close=False,
        # frameless=True,
        # easy_drag=False,
        # maximized=True
    )
    window.events.closed += on_closed
    (
        get_state,
        set_state,
        reset_state,
        new_document,
        get_document,
        find_document,
        update_document,
        remove_document,
        remove_empty_documents,
    ) = create_app_state(
        debug=os.getenv("DEBUG", "false").lower() == "true",
        platform=dict(
            platform=platform.platform(),
            cores=get_cpu_count(),
            memory=psutil.virtual_memory().total,
        ),
        on_update=lambda _: window
        and push_backend_state(window),  # Check if window is not None
    )
    webview.start(debug=get_state().debug, private_mode=False)
