import base64
import os
import shutil
import urllib.parse
from document_state import DocState
import cluster
from constants import data_file_suffixes

image_types = ["heatmap", "clustermap", "histogram", "violin"]

def find_source_files(state: DocState, prefix, suffixes):
    with os.scandir(state.tempdir_path) as entries:
        for entry in entries:
            if (
                entry.is_file()
                and entry.name.startswith(prefix)
                and any(
                    os.path.splitext(entry.name)[0].endswith(suffix)
                    for suffix in suffixes
                )
            ):
                yield entry


def save_image_from_api(data, format, destination):
    encoded = data.split(",")[1]

    if format == "svg":
        data = urllib.parse.unquote(encoded)
        with open(destination, "w", encoding="utf-8") as file:
            file.write(data)
    else:
        data = base64.b64decode(encoded)
        with open(destination, "wb") as file:
            file.write(data)

def prepare_export_data(export_path: str, matrix_path: str, doc: DocState, args: dict):
    if doc.filetype == "text/fasta":
        prefix = os.path.splitext(doc.basename)[0]
        suffixes = data_file_suffixes
    else:
        prefix = os.path.splitext(doc.basename)[0].removesuffix("_mat")
        suffixes = ["_mat"]

    if args["output_cluster"] == True:
        suffixes.append("_cluster")
        cluster.export(
            matrix_path,
            args["cluster_threshold"],
            args["cluster_method"] # Note this is not the same as the cluster_method param in the runner api
        )

    image_format = str(args["image_format"]).lower()
    saveable_formats = ["jpeg", "svg", "png"]

    if image_format not in saveable_formats:
        raise Exception(f"Expected image_format to be one of {saveable_formats}")

    base_filename = os.path.basename(os.path.splitext(doc.basename)[0]).removesuffix(
        "_fasta"
    )

    image_filenames = {
        image_type: f"{base_filename}_{image_type}.{image_format}"
        for image_type in image_types
    }

    image_destinations = {
        img_type: os.path.join(export_path, filename)
        for img_type, filename in image_filenames.items()
    }

    destination_files = [
        os.path.join(export_path, entry.name)
        for entry in find_source_files(doc, prefix, suffixes)
    ]
    destination_files.extend(image_destinations.values())

    return destination_files, image_destinations, image_format, prefix, suffixes

def do_export_data(export_path, image_destinations, image_format, doc, prefix, suffixes, args: dict):
    for entry in find_source_files(doc, prefix, suffixes):
        destination_path = os.path.join(export_path, entry.name)
        temp_destination_path = destination_path + ".tmp"
        shutil.copy2(entry.path, temp_destination_path)
        os.replace(temp_destination_path, destination_path)

    for image_type in image_types:
        save_image_from_api(
            data=args[f"{image_type}_image_data"],
            format=image_format,
            destination=image_destinations[image_type],
        )
