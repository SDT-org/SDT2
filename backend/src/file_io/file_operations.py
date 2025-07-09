import os
import json
import zipfile


def read_json_file(filename: str):
    with open(filename, "r", encoding="utf-8") as file:
        content = file.read()
        if content.strip():
            return json.loads(content)
    return None


def save_json_file(data: dict, filename: str, indent: int = 4):
    with open(filename, "w") as file:
        json.dump(data, file, indent=indent)


def pack_document(target_file, files):
    if not target_file.endswith(".sdt"):
        target_file += ".sdt"
    with zipfile.ZipFile(target_file, "w") as zf:
        for file in files:
            zf.write(
                os.path.join(os.path.dirname(target_file), file),
                os.path.basename(file),
            )


def unpack_document(input_file, output_directory):
    if not input_file.endswith(".sdt"):
        raise ValueError("Invalid file extension. Expecting '.sdt'")
    with zipfile.ZipFile(input_file, "r") as zf:
        zf.extractall(output_directory)
