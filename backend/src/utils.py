import os
import psutil
from uuid import uuid4
from hashlib import sha256
import shutil
import subprocess
import platform


def get_child_process_info():
    process = psutil.Process()
    info = []

    for child in process.children(recursive=True):
        try:
            cpu_percent = child.cpu_percent(interval=0.1)
            memory = child.memory_info().rss
            info.append((child.pid, cpu_percent, memory, child.nice()))
        except:
            pass

    return info


def make_doc_id() -> str:
    return str(uuid4())


def open_folder(path: str):
    system = platform.system().lower()

    if system == "darwin":
        cmd = ["open", path]
    elif system == "windows":
        cmd = ["explorer", path]
        shell = True
    elif system == "linux":
        for file_manager in ["xdg-open", "gnome-open", "kde-open", "gio"]:
            if shutil.which(file_manager):
                cmd = [file_manager, path]
                break
        else:
            raise OSError("Could not find a file manager to open the folder")
    else:
        raise OSError(f"Unsupported platform: {system}")

    shell = system == "windows"
    subprocess.Popen(cmd, close_fds=True, shell=shell)


def friendly_total_time(total_time):
    m, s = divmod(total_time, 60)
    return f"{int(m)}m {s:.2f}s" if m > 0 else f"{s:.2f}s"


def file_exists(path):
    return os.path.exists(os.path.join(os.path.dirname(__file__), path))
