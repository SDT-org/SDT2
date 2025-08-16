import json
import os
from config import dev_frontend_host, is_compiled


def file_exists(path):
    return os.path.exists(os.path.join(os.path.dirname(__file__), path))


def get_html_path(filename="index.html"):
    if is_compiled:
        if file_exists(f"./gui/{filename}"):
            return f"./gui/{filename}"
        raise Exception(f"{filename} not found")
    else:
        return f"{dev_frontend_host}/{filename}"


def read_json_file(filename: str):
    with open(filename, "r", encoding="utf-8") as file:
        content = file.read()
        if content.strip():
            return json.loads(content)
    return None
