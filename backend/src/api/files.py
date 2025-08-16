import mimetypes
import os
from shutil import copy
from pandas.core.frame import DataFrame
import webview
from platformdirs import user_documents_dir
from api.workflow import get_compute_stats
from app_globals import get_parsed_workflow_results
from app_settings import add_recent_file, remove_recent_file
from config import matrix_filetypes, temp_dir

from app_state import get_state, new_document, remove_empty_documents
from document_paths import build_document_paths
from export_utils import save_cols_to_csv
from save_document import unpack_document
from transformations import read_csv_matrix, similarity_triangle_to_matrix
from utils import make_doc_id
from workflow.runner import run_parse
from numpy import eye, where, nan, nanmin, nanmax


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
        if not isinstance(df, DataFrame) or df.empty:
            remove_recent_file(filepath)
            raise Exception(f"File is not a valid matrix file: {filepath}")
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


class Files:
    def open_file(self, filepath: str, doc_id: str | None = None):
        return handle_open_file(filepath=filepath, doc_id=doc_id)

    def open_file_dialog(self, doc_id: str | None = None, filepath: str | None = None):
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
        return handle_open_file(result, doc_id)

    def save_file_dialog(self, filename: str, directory: str | None = None):
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

    def select_path_dialog(self, directory: str | None = None):
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
