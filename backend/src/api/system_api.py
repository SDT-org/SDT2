import os
import sys
import json
import platform
import psutil
import webview
from platformdirs import user_documents_dir

sys.path.append(os.path.join(os.path.dirname(__file__), ".."))

from config import app_version
from app_settings import load_app_settings, save_app_settings
from utils import get_child_process_info


def close_app():
    os._exit(0)


def show_about():
    from api.webview_windows import about_window
    about_window()


def show_manual():
    from api.webview_windows import manual_window
    manual_window()


def app_config():
    return {"appVersion": app_version, "userPath": os.path.expanduser("~")}


def app_settings():
    settings = load_app_settings()
    settings["recent_files"] = [
        f for f in settings["recent_files"] if os.path.exists(f)
    ]
    export_path = settings["user_settings"].get("export_path", "")
    settings["user_settings"]["export_path"] = (
        user_documents_dir() if not os.path.exists(export_path) else export_path
    )
    save_app_settings(settings)
    return settings


def processes_info():
    return json.dumps(get_child_process_info())


def get_available_memory():
    return psutil.virtual_memory().available


def get_state(get_state_func):
    return get_state_func()._asdict()


def reset_state(reset_state_func, window, default_window_title):
    reset_state_func()
    if window:
        window.title = default_window_title
