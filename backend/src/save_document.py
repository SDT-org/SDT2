import os
import zipfile
from constants import data_file_suffixes

internal_files = ["document.json", "matrix.csv"]

def pack_document(output_file, files):
    if not output_file.endswith(".sdt"):
        output_file += ".sdt"
    with zipfile.ZipFile(output_file, "w") as zf:
        for file in files:
            basename = os.path.basename(file)

            if (
                any(
                    os.path.splitext(basename)[0].endswith(suffix)
                    for suffix in data_file_suffixes
                )
                or basename in internal_files
            ):
                zf.write(
                    os.path.join(os.path.dirname(output_file), file),
                    os.path.basename(file),
                )


def unpack_document(input_file, output_directory):
    if not input_file.endswith(".sdt"):
        raise ValueError("Invalid file extension. Expecting '.sdt'")
    with zipfile.ZipFile(input_file, "r") as zf:
        zf.extractall(output_directory)
