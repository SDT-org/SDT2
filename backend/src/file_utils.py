import json


def read_json_file(filename: str):
    with open(filename, "r", encoding="utf-8") as file:
        content = file.read()
        if content.strip():
            return json.loads(content)
    return None
