import os
import sys
import json
from typing import Dict
import urllib.parse

current_file_path = os.path.abspath(os.path.join(os.path.dirname(__file__)))
sys.path.append(os.path.join(current_file_path, "../../"))
sys.path.append(os.path.join(current_file_path, "."))

from file_utils import read_json_file
from export_utils import save_cols_to_csv
from workflow import cluster
import numpy as np
from multiprocessing import Manager, cpu_count as get_cpu_count
from tempfile import TemporaryDirectory
from platformdirs import user_documents_dir
from utils import get_child_process_info, make_doc_id, open_folder
from workflow.models import (
    LzaniSettings,
    ParasailSettings,
    RunSettings,
    WorkflowResult,
    WorkflowRun,
)
from workflow.runner import run_process, run_parse
from document_state import save_document_settings
from app_settings import (
    add_recent_file,
    load_app_settings,
    remove_recent_file,
    save_app_settings,
    update_app_settings,
)
from export import (
    ImageFormat,
    build_source_target_pairs,
    do_export,
    save_image_from_api,
)
from save_document import pack_document, unpack_document
from app_state import create_app_state
import platform
import psutil
import webview
import json
from pandas import DataFrame
from numpy import eye, where, nan, nanmin, nanmax
import mimetypes
from shutil import copy
from debug import open_doc_folder
from config import app_version, dev_frontend_host
from constants import matrix_filetypes, default_window_title
from transformations import (
    to_triangle,
    similarity_triangle_to_matrix,
    read_csv_matrix,
    read_stats_csv,
)
from document_paths import ImageKey, build_document_paths

is_compiled = "__compiled__" in globals()
is_macos = platform.system() == "Darwin"
is_windows = platform.system() == "Windows"
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
canceled = None
start_time = None  # TODO: move into workflow
workflow_runs: Dict[str, WorkflowRun] = {}
parsed_workflow_results: Dict[str, WorkflowResult] = {}


def do_cancel_run():
    global cancel_event
    if canceled:
        # It's ok if the manager has already released the resource
        try:
            canceled.value = True
        except:
            pass
    doc_id = get_state().active_run_document_id
    if doc_id:
        del workflow_runs[doc_id]
        update_document(
            doc_id,
            view="runner",
            doc_id=doc_id,
            progress=0,
            pair_progress=0,
            pair_count=0,
        )
    set_state(active_run_document_id=None)
    print("Run canceled")


def assert_window():
    assert window is not None
    return window


def get_compute_stats(workflow_result: WorkflowResult):
    max_len = workflow_result.max_sequence_length
    state = get_state()
    # TODO: this only works for parasail for now...
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
        df = read_csv_matrix(filepath)
        # Test that the file is gonna work in get_data
        df.index.tolist()
        data = df.to_numpy()
        diag_mask = eye(data.shape[0], dtype=bool)
        data_no_diag = where(diag_mask, nan, data)
        int(nanmin(data_no_diag))
        int(nanmax(data_no_diag))

        # We need a full matrix for doing things but we don't have
        # it yet because this was a .txt/.csv lower triangle matrix
        matrix_dataframe = similarity_triangle_to_matrix(df)

        # save_cols_to_csv expects distance values, not similarity
        # The matrix_dataframe now contains distance values after similarity_triangle_to_matrix
        save_cols_to_csv(matrix_dataframe, doc_paths.columns)
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
        result = run_parse(filepath)
        if result.error:
            raise Exception(result.error[0])

        parsed_workflow_results[doc_id] = result

        compute_stats = get_compute_stats(result)
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
            result_metadata={
                "is_aa": result.is_aa,
            },
        )
        if len(get_state().documents) == 1:
            remove_empty_documents()
        add_recent_file(filepath)
        return [doc_id, str(filepath)]


def get_lzani_exec_path():  # TODO: move to dynamic config? or may need to be configurable by user
    lzani_executable_name = "lz-ani.exe" if is_windows else "lz-ani"
    
    if is_compiled:
        # When compiled, the executable is in the same directory as the main exe
        # but we need to go to backend/bin relative to the exe location
        exe_dir = os.path.dirname(sys.executable)
        lzani_executable_path = os.path.join(exe_dir, "backend", "bin", lzani_executable_name)
    else:
        # When running from source
        current_script_dir = os.path.dirname(os.path.abspath(__file__))  # .../backend/src
        backend_dir = os.path.dirname(current_script_dir)  # .../backend
        bin_dir = os.path.join(backend_dir, "bin")  # .../backend/bin
        lzani_executable_path = os.path.join(bin_dir, lzani_executable_name)
    
    if not os.path.exists(lzani_executable_path):
        raise FileNotFoundError(
            f"LZ-ANI executable not found at {lzani_executable_path}"
        )
    return lzani_executable_path


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
        export_path = settings["user_settings"].get("export_path", "")
        settings["user_settings"]["export_path"] = (
            user_documents_dir() if not os.path.exists(export_path) else export_path
        )
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

    def start_workflow_run(self, args: dict):
        global canceled
        app_state = get_state()
        if app_state.debug:
            print("\nAPI args:", args)

        if app_state.active_run_document_id:
            raise Exception("Multiple runs are not supported")

        doc = get_document(args["doc_id"])

        parsed_result = parsed_workflow_results.get(doc.id)
        if not parsed_result:
            raise Exception(
                "Parsed workflow result not found. Close the document and reselect the file."
            )

        if parsed_result.error:
            raise Exception("Workflow has an error that must be resolved.")

        settings = RunSettings(
            fasta_path=doc.filename,
            doc_paths=build_document_paths(doc.tempdir_path),
            output_path=doc.tempdir_path,
            cluster_method=args.get("cluster_method", ""),
            analysis_method=args.get("analysisMethod", "parasail"),
            lzani=LzaniSettings(
                exec_path=get_lzani_exec_path(),
                score_type=args.get("lzani_score_type", "ani"),
                aw=args.get("aw"),
                am=args.get("am"),
                mal=args.get("mal"),
                msl=args.get("msl"),
                mrd=args.get("mrd"),
                mqd=args.get("mqd"),
                reg=args.get("reg"),
                ar=args.get("ar"),
            ),
            parasail=ParasailSettings(
                process_count=max(
                    min(args.get("compute_cores", 1), get_cpu_count()), 1
                ),
                scoring_matrix=args.get("scoring_matrix"),
                open_penalty=args.get("open_penalty"),
                extend_penalty=args.get("extend_penalty"),
            ),
            export_alignments=args.get("export_alignments", False),
            alignment_export_path=args.get("alignment_export_path", ""),
        )

        workflow_run = WorkflowRun(
            result=parsed_result,
            settings=settings,
            progress=0,
        )

        workflow_runs[doc.id] = workflow_run

        update_document(doc.id, view="loader")
        with Manager() as manager:
            canceled = manager.Value("b", False)
            result = run_process(workflow_run, canceled)

            if result.error:
                if result.error == "PROCESS_CANCELED":
                    # nothing to do - the state was reset in do_cancel_run
                    return
                raise Exception(f"Workflow processing step failed: {result.error}")

        set_state(active_run_document_id=None)
        del workflow_runs[doc.id]
        update_document(doc.id, view="viewer")

    def get_workflow_run_status(self, doc_id: str):
        workflow_run = workflow_runs.get(doc_id)
        if not workflow_run:
            raise Exception(f"No workflow run found for document ID: {doc_id}")
        return {"stage": workflow_run.stage, "progress": workflow_run.progress}

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
        element = assert_window().dom.get_element(selector)
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

    def save_raster_image(
        self, doc_id: str, data: str, key: ImageKey, format: ImageFormat
    ):
        doc = get_document(doc_id)
        save_image_from_api(
            doc=doc,
            data=data,
            key=key,
            format=format,
        )
        return True

    def export(self, args: dict):
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

    def get_data(self, doc_id: str):
        # Load all data in one place
        doc = get_document(doc_id)
        doc_paths = build_document_paths(doc.tempdir_path)

        # Load matrix
        df = read_csv_matrix(doc_paths.matrix)
        ids = df.index.tolist()
        update_document(doc_id, sequences_count=len(ids))
        data = df.to_numpy()

        # Load stats if available
        if os.path.exists(doc_paths.stats):
            stats_df = read_stats_csv(doc_paths.stats)
        else:
            stats_df = DataFrame([])

        # Generate identity scores from matrix for CSV imports
        # This ensures distribution plots have data to work with
        identity_scores = []
        all_scores = []
        id_map = {id: idx for idx, id in enumerate(ids)}
        
        # Check if LZANI was used (from run settings if available)
        is_lzani = False
        if file_exists(doc_paths.run_settings):
            run_settings = read_json_file(doc_paths.run_settings)
            is_lzani = run_settings.get("analysis_method") == "lzani"
        
        # The matrix contains distance values, convert to similarity for identity scores
        for i in range(len(ids)):
            for j in range(i):
                # Convert distance to similarity (identity score)
                identity_score = 100 - data[i, j]
                all_scores.append(identity_score)
                
                # For LZANI, filter out values below 70% threshold
                if not (is_lzani and identity_score < 1):
                    identity_scores.append([id_map[ids[i]], id_map[ids[j]], identity_score])
        
        # Calculate unaligned count after collecting all scores
        unaligned_count = len([s for s in all_scores if s < 3]) if is_lzani else 0
        
        # Convert matrix to triangle format for heatmap (converts distance to similarity)
        heat_data = DataFrame(data, index=ids)
        heat_data = to_triangle(heat_data)

        # Calculate min/max after conversion to similarity values
        heat_data_np = heat_data.to_numpy()
        diag_mask = eye(heat_data_np.shape[0], dtype=bool)
        heat_data_no_diag = where(diag_mask, nan, heat_data_np)
        
        # For LZ-ANI, use the lowest value above 10% (to exclude weird low values like 0.05%)
        if is_lzani:
            # Filter out values below 10% for min calculation
            valid_values = heat_data_no_diag[~np.isnan(heat_data_no_diag) & (heat_data_no_diag >= 10)]
            if valid_values.size > 0:
                min_val = int(nanmin(valid_values))
            else:
                # If no values above 10%, use the absolute minimum
                min_val = int(nanmin(heat_data_no_diag))
        else:
            min_val = int(nanmin(heat_data_no_diag))
        
        max_val = int(nanmax(heat_data_no_diag))

        parsedData = heat_data.values.tolist()

        # Load run settings if available
        metadata = dict(minVal=min_val, maxVal=max_val)
        if file_exists(doc_paths.run_settings):
            metadata["run"] = read_json_file(doc_paths.run_settings)
        
        # Add unaligned count for LZANI
        if is_lzani:
            metadata["unaligned_count"] = unaligned_count
        
        # Helper function to calculate statistics
        def calculate_stats(values):
            if len(values) == 0:
                return None
            return {
                "mean": float(np.mean(values)),
                "median": float(np.median(values)),
                "std": float(np.std(values)),
                "min": float(np.min(values)),
                "max": float(np.max(values)),
                "q1": float(np.percentile(values, 25)),
                "q3": float(np.percentile(values, 75)),
                "count": len(values)
            }
        
        # Calculate statistics for scores
        if identity_scores:
            scores_array = np.array([score[2] for score in identity_scores])
            metadata["distribution_stats"] = calculate_stats(scores_array)
        
        # Calculate statistics for GC and length if stats are available
        if not stats_df.empty:
            stats_array = np.array(stats_df.values)
            metadata["gc_stats"] = calculate_stats(stats_array[:, 1])  # GC values for all sequences in dataset
            metadata["length_stats"] = calculate_stats(stats_array[:, 2])  # Length values for all sequences in dataset

        # Return all data
        data_to_dump = dict(
            metadata=metadata,
            data=[ids] + parsedData,
            ids=ids,
            identity_scores=identity_scores,
            full_stats=stats_df.values.tolist(),
        )
        return json.dumps(data_to_dump)

    def get_clustermap_data(self, doc_id: str, threshold: float, method: str):
        doc = get_document(doc_id)

        doc_paths = build_document_paths(doc.tempdir_path)
        matrix_df = read_csv_matrix(doc_paths.matrix)
        matrix_np = matrix_df.to_numpy()
        matrix_np = np.round(matrix_np, 2)
        sorted_ids = matrix_df.index.tolist()

        # Get cluster assignments
        seqid_clusters_df = cluster.get_clusters_dataframe(
            matrix_np, method, threshold, sorted_ids
        )
        seqid_clusters_df = seqid_clusters_df.rename(
            columns={
                str(seqid_clusters_df.columns[0]): "id",
                str(seqid_clusters_df.columns[1]): "cluster",
            }
        )

        # Store original cluster numbers for color consistency
        seqid_clusters_df["id"] = seqid_clusters_df["id"].astype(str)
        seqid_clusters_df["original_cluster"] = seqid_clusters_df["cluster"]

        # Get the linkage-based order to group sequences by cluster
        new_order = cluster.get_linkage_method_order(matrix_np, method, sorted_ids, threshold)

        # Create a mapping from old index to new index
        reorder_indices = [sorted_ids.index(id) for id in new_order]

        # Reorder the matrix
        reordered_matrix = matrix_np[np.ix_(reorder_indices, reorder_indices)]

        # Reorder the tick text
        reordered_tick_text = [str(sorted_ids[i]) for i in reorder_indices]

        # Now assign sequential cluster numbers based on visual order
        # creates a mapping of sequence ID to its cluster
        id_to_cluster = {row["id"]: row["cluster"] for _, row in seqid_clusters_df.iterrows()}

        # get cluster stats  clusters
        seen_clusters = {}
        next_cluster_num = 1
        
        total_clusters= len(seqid_clusters_df["cluster"].unique())

        largest_cluster = max(seqid_clusters_df["cluster"].value_counts())
 
       
        cluster_counts = seqid_clusters_df["cluster"].value_counts()
        singleton_clusters = cluster_counts[cluster_counts == 1]
        singletons = len(singleton_clusters)
        print(singletons, "singletons")
        # Go through sequences in descending order from the topmosost cluster and assign asending sequential cluster numbers
        for seq_id in reordered_tick_text:
            original_cluster = id_to_cluster.get(seq_id)
            if original_cluster is not None and original_cluster not in seen_clusters:
                seen_clusters[original_cluster] = next_cluster_num
                next_cluster_num += 1

        # Update the dataframe with new sequential cluster numbers
        for idx, row in seqid_clusters_df.iterrows():
            original_cluster = row["cluster"]
            if original_cluster in seen_clusters:
                seqid_clusters_df.at[idx, "cluster"] = seen_clusters[original_cluster]

        # Update the cluster data to match the new order so that the legend matches
        cluster_data = seqid_clusters_df.to_dict(orient="records")

        cluster_stats = {
            "total_clusters": total_clusters,
            "largest_cluster": int(largest_cluster),
            "singleton_clusters": singletons,
        }

        return {
            "matrix": to_triangle(reordered_matrix, fill_value=None).tolist(),
            "tickText": reordered_tick_text,
            "clusterData": cluster_data,
            "cluster_stats": cluster_stats,
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
        update_document(
            args["id"],
            dataView=args["dataView"],
            heatmap=args["heatmap"],
            clustermap=args["clustermap"],
            distribution=args["distribution"],
        )
        doc = get_document(args["id"])
        save_document_settings(doc)

    def close_doc(self, doc_id: str):
        remove_document(doc_id)

    def set_window_title(self, title: str):
        assert_window().title = f"{title}{window_title_suffix}"

    def open_doc_folder(self, doc_id: str):
        doc = get_document(doc_id)
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
    if window is None:
        # This can happen if the state updates before the window is fully initialized
        print(
            "Warning: push_backend_state called with window=None. Skipping UI update."
        )
        return
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
        min_size=(900, 800),
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
        find_document,
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
        on_update=lambda _: window
        and push_backend_state(window),  # Check if window is not None
    )
    webview.start(debug=get_state().debug, private_mode=False)
