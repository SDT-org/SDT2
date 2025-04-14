import json
import os
from platformdirs import user_data_dir
current_version = 1

default_settings = {
    "recent_files": [],
    "user_settings": {}
}

def get_app_settings_path():
    data_dir = user_data_dir(appname="SDT2", ensure_exists=True)
    return os.path.join(data_dir, "settings.json")


def load_app_settings():
    path = get_app_settings_path()
    if os.path.exists(path):
        with open(path) as f:
            return json.load(f)
    else:
        return default_settings

def add_recent_file(filepath: str):
    settings = load_app_settings()
    if filepath not in settings["recent_files"]:
        settings["recent_files"].insert(0, filepath)
        settings["recent_files"] = settings["recent_files"][:8]
        save_app_settings(settings)

def remove_recent_file(filepath: str):
    settings = load_app_settings()
    if filepath in settings["recent_files"]:
        settings["recent_files"].remove(filepath)
        save_app_settings(settings)

def save_app_settings(settings: dict):
    settings["version"] = current_version
    with open(get_app_settings_path(), "w") as f:
        json.dump(settings, f, indent=2)
