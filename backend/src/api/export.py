import os
import mimetypes
import urllib.parse

from app_settings import update_app_settings
from constants import matrix_filetypes
from document_paths import ImageKey, ImageFormat, build_document_paths
from export import build_source_target_pairs, do_export, save_image_from_api
from utils import open_folder
from workflow import cluster
from app_state import get_document
from app_globals import get_app_window


class Export:
    def confirm_overwrite(self, target_files):
        files_to_overwrite = [f for f in target_files if os.path.exists(f)]
        if not files_to_overwrite:
            return True
        message = (
            "The following files will be overwritten:\n\n"
            + "\n".join(map(os.path.basename, files_to_overwrite))
            + "\n\nAre you sure?"
        )
        return get_app_window().create_confirmation_dialog("Warning", message)

    def save_svg_element(self, doc_id: str, selector: str, key: ImageKey):
        doc = get_document(doc_id)
        element = get_app_window().dom.get_element(selector)
        if element == None:
            raise Exception(f"could not find element: {selector}")
        inner_html = element.node["innerHTML"]
        data = f"<svg xmlns='http://www.w3.org/2000/svg'>{inner_html}</svg>"
        save_image_from_api(doc=doc, data=data, key=key, format="svg")
        return True

    def save_svg_data(self, doc_id: str, data: str, key: ImageKey, format: ImageFormat):
        doc = get_document(doc_id)
        data = urllib.parse.unquote(data.split(",")[1])
        save_image_from_api(
            doc=doc,
            data=data,
            key=key,
            format=format,
        )
        return True

    def save_raster_image(self, doc_id: str, data: str, key: ImageKey, format: ImageFormat):
        doc = get_document(doc_id)
        save_image_from_api(
            doc=doc,
            data=data,
            key=key,
            format=format,
        )
        return True

    def export_data(self, args: dict):
        doc = get_document(args["doc_id"])
        update_app_settings(
            {
                "user_settings": {
                    "export_path": args["export_path"],
                    "open_folder_after_export": args["open_folder"],
                }
            }
        )
        doc_paths = build_document_paths(doc.tempdir_path)
        if args["output_cluster"] == True:
            cluster.export(
                matrix_path=doc_paths.matrix,
                cluster_data_output_dir=doc_paths.cluster_dir,
                seq_dict_path=doc_paths.seq_dict,
                threshold=args["cluster_threshold"],
                method=args["cluster_method"],
            )
        prefix_default = os.path.splitext(doc.basename)[0]
        args["prefix"] = args.get("prefix", prefix_default) or prefix_default
        filetype, _ = mimetypes.guess_type(doc.basename)
        matrix_only = filetype in matrix_filetypes
        source_target_pairs = build_source_target_pairs(
            doc.tempdir_path,
            args["export_path"],
            args["prefix"],
            args["image_format"],
            matrix_only,
        )
        source_paths, target_paths = zip(*source_target_pairs)
        for file in source_paths:
            if not os.path.exists(file):
                raise FileNotFoundError(f"File not found: {file}")
        if not self.confirm_overwrite(target_paths):
            return False
        do_export(source_target_pairs)
        if args["open_folder"]:
            open_folder(args["export_path"])
        return True
