import json
import os
import psutil
from platformdirs import user_documents_dir

from app_settings import load_app_settings, save_app_settings
from config import app_version
from utils import get_child_process_info


class System:
    def get_app_config(self):
        return {"appVersion": app_version, "userPath": os.path.expanduser("~")}

    def get_app_settings(self):
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

    def get_processes_info(self):
        return json.dumps(get_child_process_info())

    def get_available_memory(self):
        return psutil.virtual_memory().available

    def close_app(self):
        os._exit(0)
