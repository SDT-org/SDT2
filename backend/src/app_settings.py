import json
import os
import time
from platformdirs import user_data_dir
current_version = 1

default_settings = {
    "recent_files": [],
    "user_settings": {
        "export_path": "",
        "open_folder_after_export": False
    }
}

def get_app_settings_path():
    data_dir = user_data_dir(appname="SDT2", ensure_exists=True)
    return os.path.join(data_dir, "settings.json")

def load_app_settings():
    path = get_app_settings_path()
    for attempt in range(3):
        if os.path.exists(path):
            try:
                with open(path, 'r', encoding='utf-8') as f:
                    content = f.read()
                    if content.strip():
                        return json.loads(content)
            except (json.JSONDecodeError, IOError) as e:
                if attempt == 2:
                    print(f"Failed to load settings after retries: {e}")

        if attempt < 2:
            time.sleep(0.05)
    else:
        return default_settings

def update_app_settings(new_settings: dict):
    settings = load_app_settings()
    settings.update(new_settings)
    save_app_settings(settings)

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
