import shutil
import subprocess
import platform

from document_state import DocState

def open_doc_folder(doc: DocState):
    path = doc.tempdir_path
    system = platform.system().lower()

    if system == 'darwin':
        cmd = ['open', path]
    elif system == 'windows':
        cmd = ['explorer', path]
        shell = True
    elif system == 'linux':
        for file_manager in ['xdg-open', 'gnome-open', 'kde-open', 'gio']:
            if shutil.which(file_manager):
                cmd = [file_manager, path]
                break
        else:
            raise OSError("Could not find a file manager to open the folder")
    else:
        raise OSError(f"Unsupported platform: {system}")

    shell = system == 'windows'
    subprocess.Popen(cmd, close_fds=True, shell=shell)
