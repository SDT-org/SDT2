import os
import sys
import json
from typing import Dict
import urllib.parse

current_file_path = os.path.abspath(os.path.join(os.path.dirname(__file__)))
sys.path.append(os.path.join(current_file_path, "../../"))
sys.path.append(os.path.join(current_file_path, "."))

from file_utils import read_json_file
from export_utils import save_cols_to_csv
from workflow import cluster
import numpy as np
from multiprocessing import Manager, cpu_count as get_cpu_count
from tempfile import TemporaryDirectory
from platformdirs import user_documents_dir
from utils import get_child_process_info, make_doc_id, open_folder
from workflow.models import (
    LzaniSettings,
    ParasailSettings,
    RunSettings,
    WorkflowResult,
    WorkflowRun,
)
from workflow.runner import run_process, run_parse
from document_state import save_document_settings
from app_settings import (
    add_recent_file,
    remove_recent_file,
)
from export import (
    ImageFormat,
)
from save_document import unpack_document
from app_state import (
    initialize_app_state,
    get_state,
    remove_empty_documents,
    set_state,
    update_document,
    new_document,
)
import app_globals
from app_globals import (
    initialize_app_globals,
    get_workflow_runs,
    get_parsed_workflow_results,
    get_canceled,
)
import platform
import psutil
import webview
import json
from numpy import eye, where, nan, nanmin, nanmax
import mimetypes
from shutil import copy
from config import dev_frontend_host
from constants import matrix_filetypes, default_window_title
from transformations import (
    similarity_triangle_to_matrix,
    read_csv_matrix,
)
from document_paths import ImageKey, build_document_paths
from mime_setup import register_mimetypes

register_mimetypes()

is_compiled = "__compiled__" in globals()
is_macos = platform.system() == "Darwin"
is_windows = platform.system() == "Windows"
temp_dir = TemporaryDirectory()
window_title_suffix = "" if is_macos else " - SDT2"
try:
    cpu_count = get_cpu_count()
except:
    cpu_count = 1
start_time = None  # TODO: move into workflow


def do_cancel_run():
    global cancel_event

    canceled = get_canceled()
    if canceled:
        # It's ok if the manager has already released the resource
        try:
            canceled.value = True
        except:
            pass
    doc_id = get_state().active_run_document_id
    if doc_id:
        workflow_runs = get_workflow_runs()
        if doc_id in workflow_runs:
            del workflow_runs[doc_id]
        update_document(
            doc_id,
            view="runner",
            doc_id=doc_id,
            progress=0,
            pair_progress=0,
            pair_count=0,
        )
    set_state(active_run_document_id=None)
    print("Run canceled")


def get_compute_stats(workflow_result: WorkflowResult):
    max_len = workflow_result.max_sequence_length
    state = get_state()
    # TODO: this only works for parasail for now...
    required_memory = (
        max_len * max_len
    ) + 100000000  # Each process has a minimum of about 100MB
    available_memory = psutil.virtual_memory().available
    total_cores = state.platform["cores"]
    min_cores = available_memory // required_memory
    if required_memory > available_memory:
        min_cores = 0
    return {
        "recommended_cores": min(
            max(round(total_cores * 0.75), 1),
            min_cores,
        ),
        "required_memory": required_memory,
        "available_memory": available_memory,
    }


def handle_open_file(filepath: str, doc_id: str | None):
    if not os.path.exists(filepath):
        remove_recent_file(filepath)
        raise Exception(f"File not found: {filepath}")
    basename = os.path.basename(filepath)
    filetype, _ = mimetypes.guess_type(basename)
    if doc_id == None:
        for doc in get_state().documents:
            if doc.filename == filepath:
                return [doc.id, doc.filename]
        doc_id = make_doc_id()
    unique_dir = os.path.join(temp_dir.name, doc_id)
    os.makedirs(unique_dir, exist_ok=True)
    doc_paths = build_document_paths(unique_dir)

    if filetype == "application/vnd.sdt":
        unpack_document(filepath, unique_dir)
        if not os.path.exists(doc_paths.matrix):
            remove_recent_file(filepath)
            raise Exception(f"File is not a valid SDT file: {filepath}")
        new_document(
            doc_id,
            view="viewer",
            filename=filepath,
            filemtime=os.path.getmtime(filepath),
            tempdir_path=unique_dir,
            basename=basename,
            pair_progress=0,
            pair_count=0,
            filetype=filetype,
        )
        add_recent_file(filepath)
        return [doc_id, filepath]

    if filetype in matrix_filetypes:
        copy(filepath, unique_dir)
        df = read_csv_matrix(filepath)
        # Test that the file is gonna work in get_data
        df.index.tolist()
        data = df.to_numpy()
        diag_mask = eye(data.shape[0], dtype=bool)
        data_no_diag = where(diag_mask, nan, data)
        int(nanmin(data_no_diag))
        int(nanmax(data_no_diag))

        # We need a full matrix for doing things but we don't have
        # it yet because this was a .txt/.csv lower triangle matrix
        matrix_dataframe = similarity_triangle_to_matrix(df)

        # save_cols_to_csv expects distance values, not similarity
        # The matrix_dataframe now contains distance values after similarity_triangle_to_matrix
        save_cols_to_csv(matrix_dataframe, doc_paths.columns)
        matrix_dataframe.to_csv(
            doc_paths.matrix,
            mode="w",
            header=False,
            index=True,
        )
        new_document(
            doc_id,
            view="viewer",
            filename=filepath,
            filemtime=os.path.getmtime(filepath),
            tempdir_path=unique_dir,
            basename=basename,
            pair_progress=0,
            pair_count=0,
            filetype=filetype,
        )
        add_recent_file(filepath)
        return [doc_id, filepath]
    else:
        result = run_parse(filepath)
        if result.error:
            raise Exception(result.error[0])

        get_parsed_workflow_results()[doc_id] = result

        compute_stats = get_compute_stats(result)
        new_document(
            doc_id,
            view="runner",
            filename=filepath,
            filetype=filetype,
            filemtime=os.path.getmtime(filepath),
            tempdir_path=unique_dir,
            basename=basename,
            validation_error_id=None,
            pair_progress=0,
            pair_count=0,
            stage="",
            progress=0,
            compute_stats=compute_stats,
            result_metadata={
                "is_aa": result.is_aa,
            },
        )
        if len(get_state().documents) == 1:
            remove_empty_documents()
        add_recent_file(filepath)
        return [doc_id, str(filepath)]


def get_lzani_exec_path():  # TODO: move to dynamic config? or may need to be configurable by user
    lzani_executable_name = "lz-ani.exe" if is_windows else "lz-ani"

    if is_compiled:
        # When compiled, the executable is in the same directory as the main exe
        # but we need to go to backend/bin relative to the exe location
        exe_dir = os.path.dirname(sys.executable)
        lzani_executable_path = os.path.join(
            exe_dir, "backend", "bin", lzani_executable_name
        )
    else:
        # When running from source
        current_script_dir = os.path.dirname(
            os.path.abspath(__file__)
        )  # .../backend/src
        backend_dir = os.path.dirname(current_script_dir)  # .../backend
        bin_dir = os.path.join(backend_dir, "bin")  # .../backend/bin
        lzani_executable_path = os.path.join(bin_dir, lzani_executable_name)

    if not os.path.exists(lzani_executable_path):
        raise FileNotFoundError(
            f"LZ-ANI executable not found at {lzani_executable_path}"
        )
    return lzani_executable_path


from api import system, files, documents, workflow, data, export, windows


class Api:
    def __init__(self):
        self.system = system.System()
        self.files = files.Files()
        self.documents = documents.Documents()
        self.workflow = workflow.Workflow()
        self.data = data.Data()
        self.export = export.Export()
        self.windows = windows.Windows()


def file_exists(path):
    return os.path.exists(os.path.join(os.path.dirname(__file__), path))


def get_html_path(filename="index.html"):
    if is_compiled:
        if file_exists(f"./gui/{filename}"):
            return f"./gui/{filename}"
        raise Exception(f"{filename} not found")
    else:
        return f"{dev_frontend_host}/{filename}"


def push_backend_state(window: webview.Window):
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
    do_cancel_run()
    os._exit(0)


def about_window():
    webview.create_window("About", get_html_path("about.html"), js_api=api)


def manual_window():
    webview.create_window("SDT2 Manual", get_html_path("manual.html"))


if __name__ == "__main__":
    api = Api()
    app_window = webview.create_window(
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
    app_window.events.closed += on_closed

    initialize_app_globals(window=app_window)

    initialize_app_state(
        debug=os.getenv("DEBUG", "false").lower() == "true",
        platform=dict(
            platform=platform.platform(),
            cores=get_cpu_count(),
            memory=psutil.virtual_memory().total,
        ),
        on_update=lambda _: push_backend_state(app_window),
    )

    webview.start(debug=get_state().debug, private_mode=False)
