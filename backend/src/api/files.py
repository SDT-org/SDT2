import os
import webview
from platformdirs import user_documents_dir


class Files:
    def open_file(self, filepath: str, doc_id: str | None = None):
        from app import handle_open_file
        return handle_open_file(filepath=filepath, doc_id=doc_id)

    def open_file_dialog(self, doc_id: str | None = None, filepath: str | None = None):
        from app import handle_open_file
        
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