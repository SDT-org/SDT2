import os
import sys
import platform
import threading
import webview
import subprocess
from subprocess import Popen, PIPE
from pathlib import Path
import re
import tempfile
import shutil
import json
import pandas as pd
import numpy as np
import multiprocessing
import base64
import urllib.parse
import shutil
import mimetypes
from app_state import create_app_state

is_nuitka = "__compiled__" in globals()

if is_nuitka:
    ext = ".exe" if platform.system() == "Windows" else ""
    bin_path = os.path.join(Path(os.path.dirname(sys.argv[0])), "bin")
    command_prefix = dict(
        sdt2=os.path.join(bin_path, "SDT2" + ext),
        cluster=os.path.join(bin_path, "cluster" + ext),
    )
else:
    python_path = sys.executable
    scripts_path = os.path.join(Path(os.path.dirname(__file__)).parents[0], "scripts")
    command_prefix = dict(
        sdt2=[sys.executable, os.path.join(scripts_path, "SDT2.py")],
        cluster=[sys.executable, os.path.join(scripts_path, "cluster.py")],
    )

temp_dir = tempfile.TemporaryDirectory()

try:
    cpu_count = multiprocessing.cpu_count()
except:
    cpu_count = 2

performance_profiles = dict(balanced=int(cpu_count / 2), best=cpu_count, low=1)
performance_profiles.setdefault("missing_key", cpu_count)

mimetypes.add_type("text/fasta", ".fasta")
mimetypes.add_type("text/fasta", ".fas")
mimetypes.add_type("text/fasta", ".faa")
mimetypes.add_type("text/fasta", ".fnt")
mimetypes.add_type("text/fasta", ".fa")
matrix_filetypes = (
    "text/csv",
    "application/vnd.ms-excel",
    "text/plain"
)

def get_matrix_path():
    state = get_state()

    if state.filetype == 'text/fasta':
        file_base = os.path.splitext(state.basename)[0]
        return os.path.join(state.tempdir_path, f"{file_base}_mat.csv")
    else:
       return state.filename[0]


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
        print(save_path)
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
        print(result)
        if not result:
            raise Exception("filename result was None")
        basename = os.path.basename(result[0])
        filetype, _ = mimetypes.guess_type(basename)
        print(basename, filetype)

        if (filetype in matrix_filetypes):
            set_state(
                view="viewer",
                filename=result,
                tempdir_path=os.path.dirname(result[0]),
                basename=basename
            )
        else:
            set_state(filename=result, filetype=filetype, basename=basename)

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

    def ls(self):
        return os.listdir(".")

    def reset_state(self):
        reset_state()

    def run_sdt2(self, args: dict):
        set_state(view="loader")
        progress_pattern = re.compile(r"progress (\d+)%")
        sequences_pattern = re.compile(r"Number of sequences: (\d+)")
        command = (
            command_prefix["sdt2"]
            if isinstance(command_prefix["sdt2"], list)
            else [command_prefix["sdt2"]]
        ) + [
            get_state().filename[0],
            "--out_dir",
            temp_dir.name,
        ]

        if args.get("alignment_type") == "local":
            command.append("--alignment_type")
            command.append("local")
        if args.get("alignment_type") == "global":
            command.append("--alignment_type")
            command.append("global")
        if args.get("cluster_method") == "Neighbor-Joining":
            command.append("--cluster_method")
            command.append("nj")
        if args.get("cluster_method") == "UPGMA":
            command.append("--cluster_method")
            command.append("upgma")
        if args.get("cluster_method") == "None":
            command.append("--cluster_method")
            command.append("None")
        if args.get("performance_profile"):
            processes = performance_profiles[str(args.get("performance_profile"))]
            command.append("--num_processes")
            command.append(str(processes))

        alignment_output_path = get_state().alignment_output_path
        match = args.get("match")
        mismatch = args.get("mismatch")
        iog = args.get("iog")
        ieg = args.get("ieg")
        log = args.get("log")
        leg = args.get("leg")
        rog = args.get("rog")
        reg = args.get("reg")

        if match:
            command.append("--match")
            command.append(str(match))
        if mismatch:
            command.append("--mismatch")
            command.append(str(mismatch))
        if iog:
            command.append("--iog")
            command.append(str(iog))
        if ieg:
            command.append("--ieg")
            command.append(str(ieg))
        if log:
            command.append("--log")
            command.append(str(log))
        if leg:
            command.append("--leg")
            command.append(str(leg))
        if rog:
            command.append("--rog")
            command.append(str(rog))
        if reg:
            command.append("--reg")
            command.append(str(reg))

        if alignment_output_path:
            command.append("--aln_out")
            command.append(str(alignment_output_path))

        if get_state().debug:
            print("\nAPI args:", args)
            print("\nRun args:", command)

        with Popen(
            command,
            stdout=PIPE,
            bufsize=1,
            universal_newlines=True,
        ) as p:
            for line in p.stdout:
                if get_state().debug:
                    print(line)

                progress_match = progress_pattern.search(line)
                sequences_match = sequences_pattern.search(line)

                if progress_match:
                    percentage = int(progress_match.group(1))
                    set_state(progress=percentage)

                if sequences_match:
                    count = int(sequences_match.group(1))
                    set_state(sequences_count=count)

            p.wait()

            if p.returncode != 0:
                raise RuntimeError(f"Subprocess exited with return code {p.returncode}")
            # TODO: return error object

            set_state(view="viewer")

    def load_data(self):
        matrix_path = get_matrix_path()

        # https://stackoverflow.com/a/57824142
        # SDT1 matrix CSVs do not have padding for columns
        with open(matrix_path, 'r') as temp_f:
            col_count = [ len(l.split(",")) for l in temp_f.readlines() ]
            column_names = [i for i in range(0, max(col_count))]

        df = pd.read_csv(matrix_path, delimiter=",", index_col=0, header=None, names=column_names)
        tickText = df.index.tolist()
        data = df.to_numpy()
        diag_mask = np.eye(data.shape[0], dtype=bool)
        data_no_diag = np.where(diag_mask, np.nan, data)
        min_val = int(np.nanmin(data_no_diag))
        max_val = int(np.nanmax(data_no_diag))
        return data, tickText, min_val, max_val

    def export_data(self, args: dict):
        state = get_state()
        matrix_path = get_matrix_path()
        suffixes = ["_cols", "_cluster", "_mat", "_summary", "_tree"]
        if state.filetype == "text/fasta":
            prefix = os.path.splitext(state.basename)[0]
        else:
            prefix = os.path.splitext(state.basename)[0].removesuffix("_mat")

        if args["output_cluster"] == True:
            command = (
                command_prefix["cluster"]
                if isinstance(command_prefix["cluster"], list)
                else [command_prefix["cluster"]]
            ) + [
                matrix_path,
                "--threshold_1",
                str(args["cluster_threshold_one"]),
                "--threshold_2",
                str(args["cluster_threshold_two"]),
            ]
            subprocess.run(command, check=True)

        with os.scandir(state.tempdir_path) as entries:
            for entry in entries:
               if entry.is_file() and entry.name.startswith(prefix) and any(os.path.splitext(entry.name)[0].endswith(suffix) for suffix in suffixes):
                   destination_path = os.path.join(state.export_path, entry.name)
                   temp_destination_path = destination_path + ".tmp"
                   shutil.copy2(entry.path, temp_destination_path)
                   os.replace(temp_destination_path, destination_path)

    def get_heatmap_data(self):
        data, tickText, min_val, max_val = self.load_data()
        heat_data = pd.DataFrame(data, index=tickText)
        parsedData = heat_data.values.tolist()
        return json.dumps([tickText] + parsedData)

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


def set_interval(interval):
    def decorator(function):
        def wrapper(*args, **kwargs):
            stopped = threading.Event()

            def loop():  # executed in another thread
                while not stopped.wait(interval):  # until stopped
                    function(*args, **kwargs)

            t = threading.Thread(target=loop)
            t.daemon = True  # stop if the program exits
            t.start()
            return stopped

        return wrapper

    return decorator


def update_client_state(window: webview.Window):
    window.evaluate_js("syncAppState()")


def on_closing():
    temp_dir.cleanup()


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
        # maximized=True
    )

    get_state, set_state, reset_state = create_app_state(
        debug=os.getenv("DEBUG", "false").lower() == "true",
        tempdir_path=temp_dir.name,
        on_update=lambda _: update_client_state(window),
    )

    webview.start(debug=get_state().debug)
