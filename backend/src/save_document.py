import os
import zipfile

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
