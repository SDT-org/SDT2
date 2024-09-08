import multiprocessing
import os
import sys
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
from app_state import create_app_state
from validations import validate_fasta
from process_data import process_data
from multiprocessing import Lock, Manager, Pool, cpu_count

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "../../")))
from config import app_version

is_nuitka = "__compiled__" in globals()
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

    set_state(
        view="runner", progress=0, pair_progress=0, pair_count=0, estimated_time=None
    )


def assert_window():
    assert window is not None
    return window


def save_image_from_api(data, format, destination):

    encoded = data.split(",")[1]

    if format == "svg":
        data = urllib.parse.unquote(encoded)
        with open(destination, "w") as file:
            file.write(data)
    else:
        data = base64.b64decode(encoded)
        with open(destination, "wb") as file:
            file.write(data)


def get_compute_stats(filename):
    max_len = 0
    for record in SeqIO.parse(filename, "fasta"):
        max_len = max(max_len, len(record.seq))

    required_memory = max_len * max_len
    total_memory = psutil.virtual_memory().total

    total_cores = multiprocessing.cpu_count()
    return {
        "total_cores": total_cores,
        "recommended_cores": min(total_cores, total_memory // required_memory),
        "total_memory": total_memory,
        "required_memory": required_memory,
        "available_memory": psutil.virtual_memory().available,
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

    def open_file_dialog(self):
        result = webview.windows[0].create_file_dialog(
            webview.OPEN_DIALOG,
            allow_multiple=False,
            file_types=(
                "Compatible file (*.fasta;*.fas;*.faa;*.fnt;*.fa;*.csv;*.txt)",
                "FASTA file (*.fasta;*.fas;*.faa;*.fnt;*.fa)",
                "SDT2 Matrix file (*.csv)",
                "SDT1 Matrix file (*.txt)",
            ),
        )
        print(result)
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
                pair_progress=0,
                pair_count=0,
                filetype=filetype,
            )
            assert_window().title = f"SDT2 - {get_state().basename}"
        else:
            assert_window().title = default_window_title
            valid, message = validate_fasta(result[0], filetype)

            if valid:
                compute_stats = get_compute_stats(result[0])
                set_state(
                    view="runner",
                    filename=result,
                    filetype=filetype,
                    basename=basename,
                    validation_error_id=None,
                    pair_progress=0,
                    pair_count=0,
                    stage="",
                    progress=0,
                    tempdir_path=temp_dir.name,
                    compute_stats=compute_stats,
                )
            else:
                set_state(
                    view="runner",
                    validation_error_id=message,
                    filename="",
                    basename="",
                    filetype="",
                    compute_stats=None,
                )

    def select_alignment_output_path(self):
        result = webview.windows[0].create_file_dialog(webview.FOLDER_DIALOG)
        print(result)
        if result:
            if isinstance(result, str):
                set_state(alignment_output_path=result)
            else:
                set_state(alignment_output_path=result[0])
        else:
            raise Exception("path result was None")

    def select_export_path(self):
        result = webview.windows[0].create_file_dialog(webview.FOLDER_DIALOG)
        if result:
            if isinstance(result, str):
                set_state(export_path=result)
            else:
                set_state(export_path=result[0])
        else:
            raise Exception("path result was None")

    def get_state(self):
        return get_state()._asdict()

    def reset_state(self):
        reset_state()
        assert_window().title = default_window_title

    def run_process_data(self, args: dict):
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

        alignment_output_path = get_state().alignment_output_path

        if alignment_output_path:
            settings["aln_out"] = str(alignment_output_path)

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
                        progress = (counter.value / pair_count) * 100
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

        if state.filetype == "text/fasta":
            prefix = os.path.splitext(state.basename)[0]
            suffixes = ["_cols", "_mat", "_summary", "_tree"]
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
        saveable_formats = ["jpg", "svg", "png"]

        if image_format not in saveable_formats:
            raise Exception(f"Expected image_format to be one of {saveable_formats}")

        heatmap_image_filename = (
            os.path.basename(os.path.splitext(state.basename)[0]).removesuffix("_fasta")
            + "_heatmap."
            + image_format
        )
        distribution_image_filename = (
            os.path.basename(os.path.splitext(state.basename)[0]).removesuffix("_fasta")
            + "_distribution."
            + image_format
        )
        heatmap_image_destination = os.path.join(
            state.export_path, heatmap_image_filename
        )
        distribution_image_destination = os.path.join(
            state.export_path, distribution_image_filename
        )

        destination_files = [
            os.path.join(state.export_path, entry.name)
            for entry in find_source_files(prefix, suffixes)
        ]
        destination_files.extend(
            [heatmap_image_destination, distribution_image_destination]
        )

        if not self.confirm_overwrite(destination_files):
            return False

        for entry in find_source_files(prefix, suffixes):
            destination_path = os.path.join(state.export_path, entry.name)
            temp_destination_path = destination_path + ".tmp"
            shutil.copy2(entry.path, temp_destination_path)
            os.replace(temp_destination_path, destination_path)

        save_image_from_api(
            data=args["heatmap_image_data"],
            format=image_format,
            destination=heatmap_image_destination,
        )
        save_image_from_api(
            data=args["distribution_image_data"],
            format=image_format,
            destination=distribution_image_destination,
        )

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


def file_exists(path):
    return os.path.exists(os.path.join(os.path.dirname(__file__), path))


def get_entrypoint(filename="index.html"):
    if file_exists(f"./gui/{filename}"):
        return f"./gui/{filename}"
    elif file_exists(f"../../gui/{filename}"):
        return f"../../gui/{filename}"

    raise Exception(f"{filename} not found")


def update_client_state(window: webview.Window):
    window.evaluate_js("syncAppState()")


def on_closed():
    temp_dir.cleanup()
    do_cancel_run()
    os._exit(0)


def about_window():
    webview.create_window("About", get_entrypoint("about.html"))


def manual_window():
    webview.create_window("SDT2 Manual", get_entrypoint("manual.html"))


entry = get_entrypoint()
if __name__ == "__main__":
    api = Api()

    window = webview.create_window(
        default_window_title,
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
    )

    # menu_items = [
    #     webview_menu.Menu(
    #         "File",
    #         [
    #             webview_menu.MenuAction("Select file...", api.open_file_dialog),
    #         ],
    #     ),
    #     webview_menu.Menu(
    #         "Help",
    #         [
    #             webview_menu.MenuAction("About SDT2", about_window),
    #             webview_menu.MenuAction("Manual", manual_window),
    #         ],
    #     ),
    # ]

    webview.start(debug=get_state().debug)
