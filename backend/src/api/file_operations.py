import os
import mimetypes
from shutil import copy
import webview
from platformdirs import user_documents_dir
from pandas import DataFrame
from numpy import eye, where, nan, nanmin, nanmax

import sys
import os as _os
sys.path.append(_os.path.join(_os.path.dirname(__file__), ".."))

from app_settings import add_recent_file, remove_recent_file
from document_paths import build_document_paths
from transformations import read_csv_matrix, similarity_triangle_to_matrix
from export_utils import save_cols_to_csv
from workflow.runner import run_parse
from constants import matrix_filetypes
from api.workflow_api import get_compute_stats


def handle_open_file(filepath: str, doc_id: str | None, temp_dir, get_state, new_document, 
                     make_doc_id, parsed_workflow_results, remove_empty_documents=None):
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
        from save_document import unpack_document
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
            raise Exception(result.error)

        parsed_workflow_results[doc_id] = result

        compute_stats = get_compute_stats(result, get_state)
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
        if remove_empty_documents and len(get_state().documents) == 1:
            remove_empty_documents()
        add_recent_file(filepath)
        return [doc_id, str(filepath)]


def open_file_dialog(doc_id: str | None, filepath: str | None, handle_open_file_func):
    if filepath is None:
        filepath = ""
    result = webview.windows[0].create_file_dialog(
        webview.OPEN_DIALOG,
        allow_multiple=False,
        directory=os.path.dirname(filepath),
        file_types=(
            "Compatible file (*.fasta;*.fas;*.faa;*.fnt;*.fa;*.sdt;*.csv;*.txt)",
            "SDT file (*.sdt)",
            "FASTA file (*.fasta;*.fas;*.faa;*.fnt;*.fa)",
            "SDT2 Matrix file (*.csv)",
            "SDT1 Matrix file (*.txt)",
        ),
    )
    if not result:
        return ""
    if isinstance(result, str):
        result = result
    else:
        result = result[0]
    return handle_open_file_func(result, doc_id)


def save_file_dialog(filename: str, directory: str | None = None):
    if directory == None:
        directory = user_documents_dir()
    result = webview.windows[0].create_file_dialog(
        webview.SAVE_DIALOG, directory=directory, save_filename=filename
    )
    output_path = None
    if result:
        if isinstance(result, str):
            output_path = result
        else:
            output_path = result[0]
    return output_path


def select_path_dialog(directory: str | None = None):
    if directory == None:
        directory = ""
    result = webview.windows[0].create_file_dialog(
        webview.FOLDER_DIALOG, directory=directory
    )
    output_path = None
    if result:
        if isinstance(result, str):
            output_path = result
        else:
            output_path = result[0]
    return output_path
