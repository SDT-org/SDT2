from utils import open_folder

from state.document_state import DocState

def open_doc_folder(doc: DocState):
    path = doc.tempdir_path
    open_folder(path)
