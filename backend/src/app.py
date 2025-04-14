import os
import sys
import time
import urllib.parse

current_file_path = os.path.abspath(os.path.join(os.path.dirname(__file__)))
sys.path.append(os.path.join(current_file_path, "../../"))
sys.path.append(os.path.join(current_file_path, "."))

from multiprocessing import Lock, Manager, Pool, cpu_count as get_cpu_count
from tempfile import TemporaryDirectory
from platformdirs import user_documents_dir
from utils import get_child_process_info, make_doc_id, open_folder
from process_data import process_data, save_cols_to_csv
from document_state import save_document_settings
from app_settings import (
    add_recent_file,
    load_app_settings,
    remove_recent_file,
    save_app_settings,
)
from export import (
    ImageFormat,
    build_source_target_pairs,
    do_export,
    save_image_from_api,
)
from save_document import pack_document, unpack_document
from app_state import create_app_state
from validations import validate_fasta
import cluster
import platform
from Bio import SeqIO
import psutil
import webview
import json
import pandas as pd
from pandas import read_csv, DataFrame
from numpy import eye, where, nan, nanmin, nanmax
import mimetypes
from time import perf_counter
from shutil import copy

from debug import open_doc_folder
from config import app_version, dev_frontend_host
from constants import matrix_filetypes, default_window_title
from heatmap import (
    dataframe_to_triangle,
    triangle_to_matrix,
    numpy_to_triangle,
)
from document_paths import ImageKey, build_document_paths

is_compiled = "__compiled__" in globals()
is_macos = platform.system() == "Darwin"
temp_dir = TemporaryDirectory()

window_title_suffix = "" if is_macos else " - SDT2"

try:
    cpu_count = get_cpu_count()
except:
    cpu_count = 1

mimetypes.add_type("text/fasta", ".fasta")
mimetypes.add_type("text/fasta", ".fas")
mimetypes.add_type("text/fasta", ".faa")
mimetypes.add_type("text/fasta", ".fnt")
mimetypes.add_type("text/fasta", ".fa")
mimetypes.add_type("application/vnd.sdt", ".sdt")


window = None
pool = None
cancelled = None
start_time = None


def do_cancel_run():
    global pool, cancelled

    if cancelled:
        # It's ok if the manager has already released the resource
        try:
            cancelled.value = True
        except:
            pass
    if pool:
        pool.close()
        pool.terminate()
        pool.join()

    set_state(active_run_document_id=None)


def assert_window():
    assert window is not None
    return window


def get_compute_stats(filename):
    max_len = 0
    for record in SeqIO.parse(filename, "fasta"):
        max_len = max(max_len, len(record.seq))

    state = get_state()

    required_memory = (
        max_len * max_len
    ) + 100000000  # Each process has a minimum of about 100MB
    available_memory = psutil.virtual_memory().available
    total_cores = state.platform["cores"]
    min_cores = available_memory // required_memory

    if required_memory > available_memory:
        min_cores = 0

    return {
        "recommended_cores": min(
            max(round(total_cores * 0.75), 1),
            min_cores,
        ),
        "required_memory": required_memory,
        "available_memory": available_memory,
    }


def handle_open_file(filepath: str, doc_id: str | None):
    if not os.path.exists(filepath):
        remove_recent_file(filepath)
        raise Exception(f"File not found: {filepath}")

    basename = os.path.basename(filepath)
    filetype, _ = mimetypes.guess_type(basename)

    if doc_id == None:
        for doc in get_state().documents:
            if doc.filename == filepath:
                return [doc.id, doc.filename]

        doc_id = make_doc_id()

    unique_dir = os.path.join(temp_dir.name, doc_id)
    os.makedirs(unique_dir, exist_ok=True)

    doc_paths = build_document_paths(unique_dir)

    if filetype == "application/vnd.sdt":
        unpack_document(filepath, unique_dir)

        if not os.path.exists(doc_paths.matrix):
            remove_recent_file(filepath)
            raise Exception(f"File is not a valid SDT file: {filepath}")

        new_document(
            doc_id,
            view="viewer",
            filename=filepath,
            filemtime=os.path.getmtime(filepath),
            tempdir_path=unique_dir,
            basename=basename,
            pair_progress=0,
            pair_count=0,
            filetype=filetype,
        )
        add_recent_file(filepath)
        return [doc_id, filepath]

    if filetype in matrix_filetypes:
        copy(filepath, unique_dir)

        # Maybe can remove this if we can find a way to make pandas
        # infer the max column length based on the latest column length
        with open(filepath, "r") as temp_f:
            col_count = [len(l.split(",")) for l in temp_f.readlines()]
            column_names = [i for i in range(0, max(col_count))]

            df = read_csv(
                filepath,
                delimiter=",",
                index_col=0,
                header=None,
                names=column_names,
            )

            # Test that the file is gonna work in get_data
            df.index.tolist()
            data = df.to_numpy()
            diag_mask = eye(data.shape[0], dtype=bool)
            data_no_diag = where(diag_mask, nan, data)
            int(nanmin(data_no_diag))
            int(nanmax(data_no_diag))

            save_cols_to_csv(df, doc_paths.triangle)

            # We need a full matrix for doing things but we don't have
            # it yet because this was a .txt/.csv lower triangle matrix
            matrix_dataframe = triangle_to_matrix(df)
            matrix_dataframe.to_csv(
                doc_paths.matrix,
                mode="w",
                header=False,
                index=True,
            )

            new_document(
                doc_id,
                view="viewer",
                filename=filepath,
                filemtime=os.path.getmtime(filepath),
                tempdir_path=unique_dir,
                basename=basename,
                pair_progress=0,
                pair_count=0,
                filetype=filetype,
            )

            add_recent_file(filepath)

            return [doc_id, filepath]

    else:
        valid, message = validate_fasta(filepath, filetype)

        if valid:
            compute_stats = get_compute_stats(filepath)
            new_document(
                doc_id,
                view="runner",
                filename=filepath,
                filetype=filetype,
                filemtime=os.path.getmtime(filepath),
                tempdir_path=unique_dir,
                basename=basename,
                validation_error_id=None,
                pair_progress=0,
                pair_count=0,
                stage="",
                progress=0,
                compute_stats=compute_stats,
            )

            if len(get_state().documents) == 1:
                remove_empty_documents()

            add_recent_file(filepath)

            return [doc_id, str(filepath)]
        else:
            raise Exception(message)


class Api:
    def close_app(self):
        os._exit(0)

    def show_about(self):
        about_window()

    def show_manual(self):
        manual_window()

    def app_config(self):
        return {"appVersion": app_version, "userPath": os.path.expanduser("~")}

    def app_settings(self):
        settings = load_app_settings()
        settings["recent_files"] = [
            f for f in settings["recent_files"] if os.path.exists(f)
        ]
        save_app_settings(settings)
        return settings

    def processes_info(self):
        return json.dumps(get_child_process_info())

    def get_available_memory(self):
        return psutil.virtual_memory().available

    def open_file(self, filepath: str, doc_id: str | None = None):
        return handle_open_file(filepath=filepath, doc_id=doc_id)

    def open_file_dialog(self, doc_id: str | None = None, filepath: str | None = None):
        if filepath is None:
            filepath = ""
        result = webview.windows[0].create_file_dialog(
            webview.OPEN_DIALOG,
            allow_multiple=False,
            directory=os.path.dirname(filepath),
            file_types=(
                "Compatible file (*.fasta;*.fas;*.faa;*.fnt;*.fa;*.sdt;*.csv;*.txt)",
                "SDT file (*.sdt)",
                "FASTA file (*.fasta;*.fas;*.faa;*.fnt;*.fa)",
                "SDT2 Matrix file (*.csv)",
                "SDT1 Matrix file (*.txt)",
            ),
        )
        if not result:
            return ""

        if isinstance(result, str):
            result = result
        else:
            result = result[0]

        return handle_open_file(result, doc_id)

    def save_file_dialog(self, filename: str, directory: str | None = None):
        if directory == None:
            directory = user_documents_dir()

        result = webview.windows[0].create_file_dialog(
            webview.SAVE_DIALOG, directory=directory, save_filename=filename
        )

        output_path = None

        if result:
            if isinstance(result, str):
                output_path = result
            else:
                output_path = result[0]

        return output_path

    def select_path_dialog(self, directory: str | None = None):
        if directory == None:
            directory = ""

        result = webview.windows[0].create_file_dialog(
            webview.FOLDER_DIALOG, directory=directory
        )
        output_path = None

        if result:
            if isinstance(result, str):
                output_path = result
            else:
                output_path = result[0]
        return output_path

    def get_state(self):
        return get_state()._asdict()

    def reset_state(self):
        reset_state()
        assert_window().title = default_window_title

    def start_run(self, args: dict):
        global pool, cancelled, start_time

        if get_state().active_run_document_id:
            raise Exception("Multiple runs are not supported")

        doc_id = args["doc_id"]
        doc = get_document(doc_id)
        if doc == None:
            raise Exception(f"could not find document: {args['doc_id']}")

        if get_state().debug:
            print("\nAPI args:", args)

        settings = dict()
        settings["input_file"] = doc.filename
        settings["out_dir"] = doc.tempdir_path

        if args.get("cluster_method") == "None":
            settings["cluster_method"] = None
        else:
            settings["cluster_method"] = args.get("cluster_method")

        compute_cores = args.get("compute_cores")
        assert isinstance(compute_cores, int)
        settings["num_processes"] = max(min(compute_cores, get_cpu_count()), 1)

        alignment_export_path = args.get("alignment_export_path")
        if args.get("export_alignments") == "True" and alignment_export_path:
            print("Saving alignments...", args.get("export_alignments"))
            settings["aln_out"] = str(alignment_export_path)

        if get_state().debug:
            print("\nRun settings:", settings)
            print(f"Number of processes: {settings['num_processes']}")

        set_state(active_run_document_id=doc_id)

        with Pool(
            settings["num_processes"],
        ) as pool:
            with Manager() as manager:
                counter = manager.Value("i", 0)
                cancelled = manager.Value("b", False)
                lock = Lock()
                start_time = perf_counter()

                def increment_pair_progress(counter, lock, start_time):
                    with lock:
                        counter.value += 1
                    doc = get_document(doc_id)
                    pair_count = doc.pair_count if doc else 0
                    if pair_count and pair_count > 0:
                        progress = round((counter.value / pair_count) * 100)
                        elapsed = perf_counter() - start_time
                        estimated_total = elapsed * (pair_count / counter.value)
                        estimated = round(estimated_total - elapsed)

                        update_document(
                            doc_id,
                            skip_callbacks=True,
                            progress=progress,
                            pair_progress=counter.value,
                            estimated_time=estimated,
                        )

                update_document(doc_id, view="loader")

                process_data(
                    settings=settings,
                    pool=pool,
                    cancelled=cancelled,
                    increment_pair_progress=lambda: increment_pair_progress(
                        counter, lock, start_time
                    ),
                    set_stage=lambda stage: update_document(doc_id, stage=stage),
                    set_pair_count=lambda pair_count: update_document(
                        doc_id, pair_count=pair_count
                    ),
                )

                if cancelled.value == False:
                    update_document(doc_id, view="viewer")
                    set_state(active_run_document_id=None)

    def cancel_run(self, doc_id: str, run_settings: str):
        do_cancel_run()

        if run_settings == "preserve":
            update_document(
                doc_id,
                view="runner",
                progress=0,
                pair_progress=0,
                pair_count=0,
                estimated_time=None,
                stage="",
                sequences_count=0,
            )
        elif run_settings == "clear":
            reset_state()

    def confirm_overwrite(self, target_files):
        files_to_overwrite = [f for f in target_files if os.path.exists(f)]

        if not files_to_overwrite:
            return True

        message = (
            "The following files will be overwritten:\n\n"
            + "\n".join(map(os.path.basename, files_to_overwrite))
            + "\n\nAre you sure?"
        )

        return assert_window().create_confirmation_dialog("Warning", message)

    def save_svg_element(self, doc_id: str, selector: str, key: ImageKey):
        doc = get_document(doc_id)
        if doc == None:
            raise Exception(f"could not find document: {doc_id}")

        element = assert_window().dom.get_element(selector)
        if element == None:
            raise Exception(f"could not find element: {selector}")
        inner_html = element.node["innerHTML"]
        data = f"<svg xmlns='http://www.w3.org/2000/svg'>{inner_html}</svg>"

        save_image_from_api(doc=doc, data=data, key=key, format="svg")

    def save_svg_data(self, doc_id: str, data: str, key: ImageKey, format: ImageFormat):
        doc = get_document(doc_id)
        if doc == None:
            raise Exception(f"could not find document: {doc_id}")

        data = urllib.parse.unquote(data.split(",")[1])
        save_image_from_api(
            doc=doc,
            data=data,
            key=key,
            format=format,
        )

    def save_raster_image(
        self, doc_id: str, data: str, key: ImageKey, format: ImageFormat
    ):
        doc = get_document(doc_id)
        if doc == None:
            raise Exception(f"could not find document: {doc_id}")

        save_image_from_api(
            doc=doc,
            data=data,
            key=key,
            format=format,
        )

    def export(self, args: dict):
        doc = get_document(args["doc_id"])
        if doc == None:
            raise Exception(f"could not find document: {args['doc_id']}")

        doc_paths = build_document_paths(doc.tempdir_path)

        if args["output_cluster"] == True:
            cluster.export(
                doc_paths.matrix,
                doc_paths.cluster,
                args["cluster_threshold"],
                # Note this is not the same as the cluster_method param in the runner api
                args["cluster_method"],
            )

        prefix_default = os.path.splitext(doc.basename)[0]
        args["prefix"] = args.get("prefix", prefix_default) or prefix_default
        source_target_pairs = build_source_target_pairs(
            doc.tempdir_path, args["export_path"], args["prefix"], args["image_format"]
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

    def load_data_and_stats(self, doc_id: str):
        doc = get_document(doc_id)
        if doc == None:
            raise Exception(f"could not find document: {doc_id}")

        doc_paths = build_document_paths(doc.tempdir_path)

        with open(doc_paths.matrix, "r") as temp_f:
            col_count = [len(l.split(",")) for l in temp_f.readlines()]
            column_names = [i for i in range(0, max(col_count))]

        if os.path.exists(doc_paths.stats):
            stats_df = read_csv(doc_paths.stats, header=0)

        else:
            stats_df = DataFrame([])

        if os.path.exists(doc_paths.columns):
            cols_data = read_csv(doc_paths.columns, skiprows=1).values.tolist()

            id_map = {}
            identity_scores = []

            for row in cols_data:
                a, b = row[:2]
                if a not in id_map:
                    id_map[a] = len(id_map)
                if b not in id_map:
                    id_map[b] = len(id_map)
                identity_scores.append([id_map[a], id_map[b]] + list(row[2:]))

            ids = list(id_map.keys())
        else:
            ids = []
            identity_scores = []

        df = read_csv(
            doc_paths.matrix,
            delimiter=",",
            index_col=0,
            header=None,
            names=column_names,
        )
        tick_text = df.index.tolist()
        update_document(doc_id, sequences_count=len(tick_text))
        data = df.to_numpy()

        diag_mask = eye(data.shape[0], dtype=bool)
        data_no_diag = where(diag_mask, nan, data)
        min_val = int(nanmin(data_no_diag))
        max_val = int(nanmax(data_no_diag))

        # TODO might be able to make one tick text object for both to use?
        return data, tick_text, min_val, max_val, ids, identity_scores, stats_df

    def get_data(self, doc_id: str):
        data, tick_text, min_val, max_val, ids, identity_scores, stats_df = (
            self.load_data_and_stats(doc_id)
        )
        heat_data = DataFrame(data, index=tick_text)
        heat_data = dataframe_to_triangle(heat_data)
        parsedData = heat_data.values.tolist()

        data_to_dump = dict(
            metadata=dict(minVal=min_val, maxVal=max_val),
            data=([tick_text] + parsedData),
            ids=ids,
            identity_scores=identity_scores,
            full_stats=stats_df.values.tolist(),
        )
        return json.dumps(data_to_dump)

    def get_clustermap_data(self, doc_id: str, threshold: float, method: str):
        doc = get_document(doc_id)
        if doc is None:
            raise Exception(f"Could not find document: {doc_id}")

        doc_paths = build_document_paths(doc.tempdir_path)

        matrix_df = pd.read_csv(
            doc_paths.matrix, delimiter=",", index_col=0, header=None
        )
        matrix_np = matrix_df.to_numpy()
        row_ids = matrix_df.index.tolist()

        # Get cluster assignments
        seqid_clusters_df = cluster.get_clusters_dataframe(
            matrix_np, method, threshold, row_ids
        )

        seqid_clusters_df = seqid_clusters_df.rename(
            columns={
                str(seqid_clusters_df.columns[0]): "id",
                str(seqid_clusters_df.columns[1]): "cluster",
            }
        )

        clusters = seqid_clusters_df["cluster"].tolist()
        reorder_clusters = cluster.order_clusters_sequentially(clusters)
        seqid_clusters_df["cluster"] = reorder_clusters

        # Sort by group
        seqid_clusters_df = seqid_clusters_df.sort_values(by=["cluster", "id"])
        sorted_ids = seqid_clusters_df["id"].tolist()

        # Get the new index order and reorder the matrix in one step using numpy
        new_order = [row_ids.index(id) for id in sorted_ids if id in row_ids]
        reordered_matrix_np = matrix_np[new_order, :][:, new_order]

        reordered_data = {
            "matrix": numpy_to_triangle(reordered_matrix_np).tolist(),
            "tickText": sorted_ids,
        }

        cluster_data = seqid_clusters_df.to_dict(orient="records")

        return {
            "matrix": reordered_data["matrix"],
            "tickText": reordered_data["tickText"],
            "clusterData": cluster_data,
        }

    def new_doc(self):
        id = make_doc_id()
        new_document(id)
        return id

    def get_doc(self, doc_id: str):
        doc = get_document(doc_id)
        return doc._asdict() if doc else None

    def save_doc(self, doc_id: str, path: str, save_as: bool = False):
        doc = get_document(doc_id)
        if doc == None:
            raise Exception(f"Expected to find document: {doc_id}")

        files = [
            os.path.join(doc.tempdir_path, file)
            for file in os.listdir(doc.tempdir_path)
        ]

        pack_document(path, files)

        if save_as == False:
            basename = os.path.basename(path)
            update_document(
                doc_id, filename=path, basename=basename, filetype="application/vnd.sdt"
            )
            add_recent_file(path)

    def save_doc_settings(self, args: dict):
        doc = get_document(args["id"])
        if doc == None:
            raise Exception(f"Expected to find document: {args['id']}")
        update_document(
            args["id"],
            dataView=args["dataView"],
            heatmap=args["heatmap"],
            clustermap=args["clustermap"],
            distribution=args["distribution"],
        )
        doc = get_document(args["id"])
        if doc == None:
            raise Exception(f"Expected to find document: {args['id']}")
        save_document_settings(doc)

    def close_doc(self, doc_id: str):
        remove_document(doc_id)

    def set_window_title(self, title: str):
        assert_window().title = f"{title}{window_title_suffix}"

    def open_doc_folder(self, doc_id: str):
        doc = get_document(doc_id)
        if doc == None:
            raise Exception(f"Expected to find document: {doc_id}")
        return open_doc_folder(doc)


def file_exists(path):
    return os.path.exists(os.path.join(os.path.dirname(__file__), path))


def get_html_path(filename="index.html"):
    if is_compiled:
        if file_exists(f"./gui/{filename}"):
            return f"./gui/{filename}"
        raise Exception(f"{filename} not found")
    else:
        return f"{dev_frontend_host}/{filename}"


def push_backend_state(window: webview.Window):
    state = get_state()
    dict_state = lambda t: {
        f: [d._asdict() for d in getattr(t, f)] if f == "documents" else getattr(t, f)
        for f in t._fields
    }
    js_app_state = json.dumps(dict(state=dict_state(state)))
    window.evaluate_js(
        f"document.dispatchEvent(new CustomEvent('sync-state', {{ detail: {js_app_state} }}))"
    )


def on_closed():
    # map(lambda doc: doc.cleanup(), get_state().documents)
    do_cancel_run()
    os._exit(0)


def about_window():
    webview.create_window("About", get_html_path("about.html"), js_api=api)


def manual_window():
    webview.create_window("SDT2 Manual", get_html_path("manual.html"))


if __name__ == "__main__":
    api = Api()

    window = webview.create_window(
        default_window_title,
        url=get_html_path(),
        js_api=api,
        # TODO: store last window size and position
        width=1200,
        height=900,
        min_size=(640, 480),
        confirm_close=False,
        # frameless=True,
        # easy_drag=False,
        # maximized=True
    )

    window.events.closed += on_closed

    (
        get_state,
        set_state,
        reset_state,
        new_document,
        get_document,
        update_document,
        remove_document,
        remove_empty_documents,
    ) = create_app_state(
        debug=os.getenv("DEBUG", "false").lower() == "true",
        platform=dict(
            platform=platform.platform(),
            cores=get_cpu_count(),
            memory=psutil.virtual_memory().total,
        ),
        on_update=lambda _: push_backend_state(window),
    )

    webview.start(debug=get_state().debug, private_mode=False)
