import os
import webview
import subprocess
import tempfile
import shutil
import json
import pandas as pd
import numpy as np
import base64
import urllib.parse
import shutil
import mimetypes
import math
import cluster
from warnings import warn
from time import perf_counter
from app_state import create_app_state
from validations import validate_fasta
from process_data import process_data
from multiprocessing import Lock, Manager, Pool, cpu_count

is_nuitka = "__compiled__" in globals()
window = None
temp_dir = tempfile.TemporaryDirectory()

try:
    cpu_count = cpu_count()
except:
    cpu_count = 1

performance_profiles = dict(
    balanced=int(math.floor(max([cpu_count / 2, 1]))),
    best=cpu_count,
    high=max(cpu_count - 1, 1),
    low=1,
)
performance_profiles.setdefault("missing_key", cpu_count)

mimetypes.add_type("text/fasta", ".fasta")
mimetypes.add_type("text/fasta", ".fas")
mimetypes.add_type("text/fasta", ".faa")
mimetypes.add_type("text/fasta", ".fnt")
mimetypes.add_type("text/fasta", ".fa")
matrix_filetypes = ("text/csv", "application/vnd.ms-excel", "text/plain")

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
    else:
        warn("Expected pool instance")

    set_state(
        view="runner", progress=0, pair_progress=0, pair_count=0, estimated_time=None
    )


class Api:
    def fullscreen(self):
        webview.windows[0].toggle_fullscreen()

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

        encoded = args["data"].split(",")[1]

        if args["format"] == "svg":
            data = urllib.parse.unquote(encoded)
            with open(save_path, "w") as file:
                file.write(data)
        else:
            data = base64.b64decode(encoded)
            with open(save_path, "wb") as file:
                file.write(data)

    def save_file_dialog(self):
        state = get_state()
        save_filename = f"{os.path.splitext(state.basename)[0]}.svg"
        save_path = webview.windows[0].create_file_dialog(
            dialog_type=webview.SAVE_DIALOG,
            directory=state.filename[0],
            save_filename=save_filename,
        )

        if not save_path:
            return

        basename = state.basename
        image_path = os.path.join(state.tempdir_path, f"{basename}.svg")

        try:
            shutil.copy2(image_path, save_path)
            print(f"File successfully copied to {save_path}")
        except IOError as e:
            print(f"An error occurred while copying the file: {e}")

    def open_file_dialog(self):
        result = webview.windows[0].create_file_dialog(
            webview.OPEN_DIALOG, allow_multiple=False
        )

        if not result:
            raise Exception("filename result was None")
        basename = os.path.basename(result[0])
        filetype, _ = mimetypes.guess_type(basename)

        if filetype in matrix_filetypes:
            set_state(
                view="viewer",
                filename=result,
                tempdir_path=os.path.dirname(result[0]),
                basename=basename,
            )
        else:
            valid, message = validate_fasta(result[0])

            if valid:
                set_state(
                    filename=result,
                    filetype=filetype,
                    basename=basename,
                    validation_error_id=None,
                )
            else:
                set_state(
                    validation_error_id=message, filename="", basename="", filetype=""
                )

    def select_alignment_output_path(self):
        result = webview.windows[0].create_file_dialog(webview.FOLDER_DIALOG)
        if result:
            set_state(alignment_output_path=result[0])
        else:
            raise Exception("path result was None")

    def select_export_path(self):
        result = webview.windows[0].create_file_dialog(webview.FOLDER_DIALOG)
        if result:
            set_state(export_path=result[0])
        else:
            raise Exception("path result was None")

    def get_state(self):
        return get_state()._asdict()

    def reset_state(self):
        reset_state()

    def run_process_data(self, args: dict):
        global pool, cancelled, start_time
        set_state(view="loader")

        settings = dict()
        settings["input_file"] = get_state().filename[0]
        settings["out_dir"] = temp_dir.name

        if args.get("alignment_type") == "local":
            settings["alignment_type"] = "local"
        if args.get("alignment_type") == "global":
            settings["alignment_type"] = "global"
        if args.get("cluster_method") == "Neighbor-Joining":
            settings["cluster_method"] = "nj"
        if args.get("cluster_method") == "UPGMA":
            settings["cluster_method"] = "upgma"
        if args.get("cluster_method") == "None":
            settings["cluster_method"] = "none"
        if args.get("performance_profile"):
            settings["num_processes"] = performance_profiles[
                str(args.get("performance_profile"))
            ]

        alignment_output_path = get_state().alignment_output_path

        params = ["match", "mismatch", "iog", "ieg", "log", "leg", "rog", "reg"]
        for param in params:
            value = args.get(param)
            if value is not None:
                settings[param] = int(value)

        if alignment_output_path:
            settings["aln_out"] = str(alignment_output_path)

        if get_state().debug:
            print("\nAPI args:", args)
            print("\nRun settings:", settings)

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
                        progress = (counter.value / pair_count) * 100
                        elapsed = perf_counter() - start_time
                        estimated_total = elapsed * (pair_count / counter.value)
                        estimated = round(estimated_total - elapsed)
                        set_state(
                            progress=progress,
                            pair_progress=counter.value,
                            estimated_time=estimated,
                        )

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

    def cancel_run(self):
        do_cancel_run()

    def load_data(self):
        matrix_path = get_matrix_path()

        # https://stackoverflow.com/a/57824142
        # SDT1 matrix CSVs do not have padding for columns
        with open(matrix_path, "r") as temp_f:
            col_count = [len(l.split(",")) for l in temp_f.readlines()]
            column_names = [i for i in range(0, max(col_count))]

        df = pd.read_csv(
            matrix_path, delimiter=",", index_col=0, header=None, names=column_names
        )
        tickText = df.index.tolist()
        count = len(tickText)
        set_state(sequences_count=count)
        data = df.to_numpy()
        diag_mask = np.eye(data.shape[0], dtype=bool)
        data_no_diag = np.where(diag_mask, np.nan, data)
        min_val = int(np.nanmin(data_no_diag))
        max_val = int(np.nanmax(data_no_diag))
        return data, tickText, min_val, max_val

    def confirm_overwrite(self, destination_files):
        if not window:
            return

        files_to_overwrite = [f for f in destination_files if os.path.exists(f)]

        if not files_to_overwrite:
            return True

        message = (
            "The following files will be overwritten:\n\n"
            + "\n".join(map(os.path.basename, files_to_overwrite))
            + "\n\nAre you sure?"
        )

        return window.create_confirmation_dialog("Warning", message)

    def export_data(self, args: dict):
        state = get_state()
        matrix_path = get_matrix_path()

        if state.filetype == "text/fasta":
            prefix = os.path.splitext(state.basename)[0]
            suffixes = ["_cols", "_mat", "_summary", "_tree"]
        else:
            prefix = os.path.splitext(state.basename)[0].removesuffix("_mat")
            suffixes = ["_mat"]

        if args["output_cluster"] == True:
            suffixes.append("_cluster")

        if args["output_cluster"] == True:
            cluster.export(
                matrix_path,
                args["cluster_threshold_one"],
                args["cluster_threshold_two"],
            )

        destination_files = [
            os.path.join(state.export_path, entry.name)
            for entry in find_source_files(prefix, suffixes)
        ]

        if not self.confirm_overwrite(destination_files):
            return False

        for entry in find_source_files(prefix, suffixes):
            destination_path = os.path.join(state.export_path, entry.name)
            temp_destination_path = destination_path + ".tmp"
            shutil.copy2(entry.path, temp_destination_path)
            os.replace(temp_destination_path, destination_path)

        return True

    def get_heatmap_data(self):
        data, tickText, min_val, max_val = self.load_data()
        heat_data = pd.DataFrame(data, index=tickText)
        parsedData = heat_data.values.tolist()
        return json.dumps(
            dict(
                metadata=dict(minVal=min_val, maxVal=max_val),
                data=([tickText] + parsedData),
            )
        )

    def get_line_histo_data(self):
        # caluclating hist data manually to allow for line and scatter
        data, _, min_val, max_val = self.load_data()
        # Create a mask for the diagonal elements
        diag_mask = np.eye(data.shape[0], dtype=bool)

        # Exclude the diagonal elements
        data_no_diag = np.where(diag_mask, np.nan, data)

        # Flatten the data and remove NaN values
        clean_vals = data_no_diag[~np.isnan(data_no_diag)].flatten()
        clean_vals = np.rint(clean_vals)
        # get unique values count number of occurances
        uniquev, countsv = np.unique(clean_vals, return_counts=True)
        props = np.around(countsv / countsv.sum(), 2)
        # subtract 1 from minimum value for chart clarity
        min_val = min_val - 1
        max_val = max_val + 2

        # create list of range values with 102 as max to ensure the 100 values fully display
        range_values = np.arange(min_val, max_val, 1)

        # create empty dict of with numerical keys populated from range_values
        proportion_dict = {int(val): float(0) for val in range_values}
        # populate the values of the dictionary  with proportions stores in props as the correstpond to uniquev
        proportion_dict.update(
            {int(val): float(prop) for val, prop in zip(uniquev, props)}
        )

        data_to_dump = {
            "x": list(proportion_dict.keys()),
            "y": list(proportion_dict.values()),
        }

        return json.dumps(data_to_dump)


def get_entrypoint():
    def exists(path):
        return os.path.exists(os.path.join(os.path.dirname(__file__), path))

    if exists("../../gui/index.html"):  # unfrozen development
        return "../../gui/index.html"

    if exists("../../Resources/gui/index.html"):  # frozen py2app
        return "../../Resources/gui/index.html"

    if exists("./gui/index.html"):
        return "./gui/index.html"

    raise Exception("No index.html found")


def update_client_state(window: webview.Window):
    window.evaluate_js("syncAppState()")


def on_closed():
    temp_dir.cleanup()
    do_cancel_run()


entry = get_entrypoint()

if __name__ == "__main__":
    api = Api()
    window = webview.create_window(
        "Sequence Demarcation Tool 2 Beta",
        entry,
        js_api=api,
        # TODO: store last window size and position
        width=1200,
        height=900,
        min_size=(640, 480),
        confirm_close=True,
        # maximized=True
    )

    window.events.closed += on_closed

    get_state, set_state, reset_state = create_app_state(
        debug=os.getenv("DEBUG", "false").lower() == "true",
        tempdir_path=temp_dir.name,
        on_update=lambda _: update_client_state(window),
        performance_profiles=performance_profiles,
    )

    webview.start(debug=get_state().debug)
