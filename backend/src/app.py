import multiprocessing
import os
import sys

current_file_path = os.path.abspath(os.path.join(os.path.dirname(__file__)))
sys.path.append(os.path.join(current_file_path, "../../"))
sys.path.append(os.path.join(current_file_path, "."))

import platform
from Bio import SeqIO
import psutil
import webview
import tempfile
import shutil
import json
import pandas as pd
import numpy as np
import base64
import urllib.parse
import shutil
import mimetypes
import cluster
from time import perf_counter
from multiprocessing import Lock, Manager, Pool, cpu_count

from app_state import create_app_state
from validations import validate_fasta
from config import app_version
from utils import get_child_process_info
from process_data import process_data, save_cols_to_csv

dev_frontend_host = "http://localhost:5173"
is_compiled = "__compiled__" in globals()
is_macos = platform.system() == "Darwin"
temp_dir = tempfile.TemporaryDirectory()

default_window_title = "Sequence Demarcation Tool 2 Beta"

try:
    cpu_count = cpu_count()
except:
    cpu_count = 1

mimetypes.add_type("text/fasta", ".fasta")
mimetypes.add_type("text/fasta", ".fas")
mimetypes.add_type("text/fasta", ".faa")
mimetypes.add_type("text/fasta", ".fnt")
mimetypes.add_type("text/fasta", ".fa")
matrix_filetypes = ("text/csv", "application/vnd.ms-excel", "text/plain")

window = None
pool = None
cancelled = None
start_time = None


def get_matrix_path():
    state = get_state()

    if state.filetype == "text/fasta":
        file_base = os.path.splitext(state.basename)[0]
        return os.path.join(state.tempdir_path, f"{file_base}_mat.csv")
    else:
        return state.filename[0]


def find_source_files(prefix, suffixes):
    state = get_state()
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


def assert_window():
    assert window is not None
    return window


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


class Api:
    def close_app(self):
        os._exit(0)

    def show_about(self):
        about_window()

    def show_manual(self):
        manual_window()

    def app_config(self):
        return json.dumps({"appVersion": app_version})

    def processes_info(self):
        return json.dumps(get_child_process_info())

    def get_available_memory(self):
        return psutil.virtual_memory().available

    def save_image(self, args: dict):
        state = get_state()
        file_name = f"{state.basename}.{args['format']}"

        save_path = webview.windows[0].create_file_dialog(
            dialog_type=webview.SAVE_DIALOG,
            directory=state.filename[0],
            save_filename=file_name,
        )

        if not save_path:
            return

        save_image_from_api(
            data=args["data"],
            format=args["format"],
            destination=save_path,
        )

    def open_file_dialog(self, filepath: str | None = None):
        if filepath is None:
            filepath = ""
        result = webview.windows[0].create_file_dialog(
            webview.OPEN_DIALOG,
            allow_multiple=False,
            directory=os.path.dirname(filepath),
            file_types=(
                "Compatible file (*.fasta;*.fas;*.faa;*.fnt;*.fa;*.csv;*.txt)",
                "FASTA file (*.fasta;*.fas;*.faa;*.fnt;*.fa)",
                "SDT2 Matrix file (*.csv)",
                "SDT1 Matrix file (*.txt)",
            ),
        )
        if not result:
            return ""

        if isinstance(result, str):
            result = [result]

        matrix_path = result[0]
        basename = os.path.basename(matrix_path)

        filetype, _ = mimetypes.guess_type(basename)

        do_cancel_run()

        if filetype in matrix_filetypes:
            # Maybe can remove this if we can find a way to make pandas
            # infer the max column length based on the latest column length
            with open(matrix_path, "r") as temp_f:
                col_count = [len(l.split(",")) for l in temp_f.readlines()]
                column_names = [i for i in range(0, max(col_count))]

            df = pd.read_csv(
                matrix_path,
                delimiter=",",
                index_col=0,
                header=None,
                names=column_names,
            )

            save_cols_to_csv(
                df,
                os.path.join(
                    temp_dir.name,
                    os.path.splitext(basename)[0].removesuffix("_mat"),
                ),
            )

            set_state(
                view="viewer",
                filename=result,
                filemtime=os.path.getmtime(result[0]),
                tempdir_path=os.path.dirname(result[0]),
                basename=basename,
                pair_progress=0,
                pair_count=0,
                filetype=filetype,
            )
            assert_window().title = f"SDT2 - {get_state().basename}"

        else:
            assert_window().title = default_window_title
            valid, message = validate_fasta(matrix_path, filetype)

            if valid:
                compute_stats = get_compute_stats(result[0])
                set_state(
                    view="runner",
                    filename=result,
                    filetype=filetype,
                    filemtime=os.path.getmtime(matrix_path),
                    basename=basename,
                    validation_error_id=None,
                    pair_progress=0,
                    pair_count=0,
                    stage="",
                    progress=0,
                    tempdir_path=temp_dir.name,
                    compute_stats=compute_stats,
                )
                return str(matrix_path)
            else:
                set_state(
                    view="runner",
                    validation_error_id=message,
                    filename="",
                    basename="",
                    filetype="",
                    filemtime=None,
                    compute_stats=None,
                )

    def select_path_dialog(self, directory: str | None = None):
        if directory is None:
            directory = ""
        print("select_path_dialog", directory)
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
        set_state(view="loader")

        if get_state().debug:
            print("\nAPI args:", args)

        settings = dict()
        settings["input_file"] = get_state().filename[0]
        settings["out_dir"] = temp_dir.name

        if args.get("cluster_method") == "Neighbor-Joining":
            settings["cluster_method"] = "nj"
        if args.get("cluster_method") == "UPGMA":
            settings["cluster_method"] = "upgma"
        if args.get("cluster_method") == "None":
            settings["cluster_method"] = "none"

        compute_cores = args.get("compute_cores")
        assert isinstance(compute_cores, int)
        settings["num_processes"] = max(
            min(compute_cores, multiprocessing.cpu_count()), 1
        )

        alignment_export_path = args.get("alignment_export_path")
        if args.get("export_alignments") == "True" and alignment_export_path:
            print("Saving alignments...", args.get("export_alignments"))
            settings["aln_out"] = str(alignment_export_path)

        if get_state().debug:
            print("\nRun settings:", settings)
            print(f"Number of processes: {settings['num_processes']}")

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
                    pair_count = get_state().pair_count
                    if pair_count and pair_count > 0:
                        progress = round((counter.value / pair_count) * 100)
                        elapsed = perf_counter() - start_time
                        estimated_total = elapsed * (pair_count / counter.value)
                        estimated = round(estimated_total - elapsed)
                        set_state(
                            progress=progress,
                            pair_progress=counter.value,
                            estimated_time=estimated,
                        )

                assert_window().title = f"SDT2 - Analyzing {get_state().basename}"
                process_data(
                    settings=settings,
                    pool=pool,
                    cancelled=cancelled,
                    increment_pair_progress=lambda: increment_pair_progress(
                        counter, lock, start_time
                    ),
                    set_stage=lambda stage: set_state(stage=stage),
                    set_pair_count=lambda pair_count: set_state(pair_count=pair_count),
                )

                if cancelled.value == False:
                    set_state(view="viewer")
                    assert_window().title = f"SDT2 - {get_state().basename}"

    def cancel_run(self, run_settings: str):
        do_cancel_run()

        if run_settings == "preserve":
            set_state(
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

    def confirm_overwrite(self, destination_files):
        files_to_overwrite = [f for f in destination_files if os.path.exists(f)]

        if not files_to_overwrite:
            return True

        message = (
            "The following files will be overwritten:\n\n"
            + "\n".join(map(os.path.basename, files_to_overwrite))
            + "\n\nAre you sure?"
        )

        return assert_window().create_confirmation_dialog("Warning", message)

    def export_data(self, args: dict):
        state = get_state()
        matrix_path = get_matrix_path()
        export_path = args["export_path"]

        if state.filetype == "text/fasta":
            prefix = os.path.splitext(state.basename)[0]
            suffixes = ["_cols", "_mat", "_summary", "_tree", "_stats"]
        else:
            prefix = os.path.splitext(state.basename)[0].removesuffix("_mat")
            suffixes = ["_mat"]

        if args["output_cluster"] == True:
            suffixes.append("_cluster")
            cluster.export(
                matrix_path,
                args["cluster_threshold_one"],
                args["cluster_threshold_two"],
            )

        image_format = str(args["image_format"]).lower()
        saveable_formats = ["jpeg", "svg", "png"]

        if image_format not in saveable_formats:
            raise Exception(f"Expected image_format to be one of {saveable_formats}")

        base_filename = os.path.basename(
            os.path.splitext(state.basename)[0]
        ).removesuffix("_fasta")
        image_types = ["heatmap", "histogram", "violin", "raincloud"]
        image_filenames = {
            img_type: f"{base_filename}_{img_type}.{image_format}"
            for img_type in image_types
        }

        image_destinations = {
            img_type: os.path.join(export_path, filename)
            for img_type, filename in image_filenames.items()
        }

        destination_files = [
            os.path.join(export_path, entry.name)
            for entry in find_source_files(prefix, suffixes)
        ]
        destination_files.extend(image_destinations.values())

        if not self.confirm_overwrite(destination_files):
            return False

        for entry in find_source_files(prefix, suffixes):
            destination_path = os.path.join(export_path, entry.name)
            temp_destination_path = destination_path + ".tmp"
            shutil.copy2(entry.path, temp_destination_path)
            os.replace(temp_destination_path, destination_path)

        for img_type in image_types:
            save_image_from_api(
                data=args[f"{img_type}_image_data"],
                format=image_format,
                destination=image_destinations[img_type],
            )

        return True

    def load_data_and_stats(self):
        state = get_state()
        file_base = os.path.splitext(state.basename)[0]
        matrix_path = get_matrix_path()

        with open(matrix_path, "r") as temp_f:
            col_count = [len(l.split(",")) for l in temp_f.readlines()]
            column_names = [i for i in range(0, max(col_count))]

        if state.filetype == "text/fasta":
            stats_path = os.path.join(state.tempdir_path, f"{file_base}_stats.csv")
            stats_df = pd.read_csv(stats_path, header=0)
            gc_stats = stats_df["GC %"].map(lambda value: round(value * 100)).tolist()
            len_stats = stats_df["Sequence Length"].tolist()
        elif state.filetype in matrix_filetypes:
            gc_stats = []
            len_stats = []
        else:
            raise Exception("unsupported file type")

        cols_dir = (
            state.tempdir_path if state.filetype == "text/fasta" else temp_dir.name
        )
        cols_file_base = file_base.removesuffix("_mat")
        cols_path = os.path.join(cols_dir, f"{cols_file_base}_cols.csv")
        identity_scores = pd.read_csv(cols_path, skiprows=1).values.tolist()

        df = pd.read_csv(
            matrix_path, delimiter=",", index_col=0, header=None, names=column_names
        )
        tick_text = df.index.tolist()
        set_state(sequences_count=len(tick_text))
        data = df.to_numpy()

        diag_mask = np.eye(data.shape[0], dtype=bool)
        data_no_diag = np.where(diag_mask, np.nan, data)
        min_val = int(np.nanmin(data_no_diag))
        max_val = int(np.nanmax(data_no_diag))

        # TODO might be able to make one tick text object for both to use?
        return data, tick_text, min_val, max_val, gc_stats, len_stats, identity_scores

    def get_data(self):
        data, tick_text, min_val, max_val, gc_stats, len_stats, identity_scores = (
            self.load_data_and_stats()
        )
        heat_data = pd.DataFrame(data, index=tick_text)
        parsedData = heat_data.values.tolist()

        data_to_dump = dict(
            metadata=dict(minVal=min_val, maxVal=max_val),
            data=([tick_text] + parsedData),
            identity_scores=identity_scores,
            gc_stats=list(gc_stats),
            length_stats=list(len_stats),
        )
        return json.dumps(data_to_dump)


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
    js_app_state = json.dumps(dict(state=get_state()._asdict()))
    window.evaluate_js(
        f"document.dispatchEvent(new CustomEvent('sync-state', {{ detail: {js_app_state} }}))"
    )


def on_closed():
    temp_dir.cleanup()
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

    get_state, set_state, reset_state = create_app_state(
        debug=os.getenv("DEBUG", "false").lower() == "true",
        tempdir_path=temp_dir.name,
        platform=dict(
            platform=platform.platform(),
            cores=multiprocessing.cpu_count(),
            memory=psutil.virtual_memory().total,
        ),
        on_update=lambda _: push_backend_state(window),
    )

    webview.start(debug=get_state().debug, private_mode=False)
