import os
import sys

current_file_path = os.path.abspath(os.path.join(os.path.dirname(__file__)))
sys.path.append(os.path.join(current_file_path, "../../"))
sys.path.append(os.path.join(current_file_path, "."))

import multiprocessing
from tempfile import TemporaryDirectory
from platformdirs import user_documents_dir
from utils import get_child_process_info, make_doc_id
from process_data import process_data, save_cols_to_csv
from document_state import DocState, save_doc_settings
from app_settings import add_recent_file, load_app_settings
from export_data import do_export_data, prepare_export_data
from save_document import pack_document, unpack_document
from app_state import create_app_state
from validations import validate_fasta
import platform
from Bio import SeqIO
import psutil
import webview
import json
from pandas import read_csv, DataFrame
from numpy import eye, where, nan, nanmin, nanmax
import mimetypes
from time import perf_counter
from multiprocessing import Lock, Manager, Pool, cpu_count

from config import app_version, dev_frontend_host
from constants import matrix_filetypes, default_window_title

is_compiled = "__compiled__" in globals()
is_macos = platform.system() == "Darwin"
temp_dir = TemporaryDirectory()

window_title_suffix = "" if is_macos else " - SDT2"

try:
    cpu_count = cpu_count()
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


def get_matrix_path(state: DocState):
    if state.filetype == "application/vnd.sdt":
        path = next(
            (
                os.path.join(state.tempdir_path, file)
                for file in os.listdir(state.tempdir_path)
                if file.endswith("_mat.csv")
            ),
            None,
        )
        if path == None:
            raise Exception(f"Failed to located matrix file")
        return path
    elif state.filetype == "text/fasta":
        file_base = os.path.splitext(state.basename)[0]
        return os.path.join(state.tempdir_path, f"{file_base}_mat.csv")
    else:
        return state.filename


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
    basename = os.path.basename(filepath)
    filetype, _ = mimetypes.guess_type(basename)

    if doc_id == None:
        for doc in get_state().documents:
            if doc.filename == filepath:
                return [doc.id, doc.filename]

        doc_id = make_doc_id()

    unique_dir = os.path.join(temp_dir.name, doc_id)
    os.makedirs(unique_dir, exist_ok=True)

    if filetype == "application/vnd.sdt":
        unpack_document(filepath, unique_dir)
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

        save_cols_to_csv(
            df,
            os.path.join(
                unique_dir,
                os.path.splitext(basename)[0].removesuffix("_mat"),
            ),
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
        return json.dumps({"appVersion": app_version})

    def app_settings(self):
        return load_app_settings()

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
        doc = get_document(args["doc_id"])
        if doc == None:
            raise Exception(f"could not find document: {args['doc_id']}")

        matrix_path = get_matrix_path(doc)
        export_path = args["export_path"]
        destination_files, image_destinations, image_format, prefix, suffixes = (
            prepare_export_data(export_path, matrix_path, doc, args)
        )

        if not self.confirm_overwrite(destination_files):
            return False

        do_export_data(
            export_path, image_destinations, image_format, doc, prefix, suffixes, args
        )

        return True

    def load_data_and_stats(self, doc_id: str):
        doc = get_document(doc_id)
        if doc == None:
            raise Exception(f"could not find document: {doc_id}")

        if doc.filetype == "application/vnd.sdt" or doc.filetype == "text/fasta":
            matrix_path = get_matrix_path(doc)
        else:
            matrix_path = doc.filename

        if doc.filetype == "application/vnd.sdt":
            file_base = os.path.splitext(os.path.basename(matrix_path))[0].removesuffix(
                "_mat"
            )
        else:
            file_base = os.path.splitext(doc.basename)[0]

        with open(matrix_path, "r") as temp_f:
            col_count = [len(l.split(",")) for l in temp_f.readlines()]
            column_names = [i for i in range(0, max(col_count))]

        if doc.filetype == "text/fasta" or doc.filetype == "application/vnd.sdt":
            stats_path = os.path.join(doc.tempdir_path, f"{file_base}_stats.csv")
            stats_df = read_csv(stats_path, header=0)

        elif doc.filetype in matrix_filetypes:
            stats_df = DataFrame([])
        else:
            raise Exception("unsupported file type")

        # cols_dir = doc.tempdir_path if doc.filetype == "text/fasta" else temp_dir.name
        cols_dir = doc.tempdir_path
        cols_file_base = file_base.removesuffix("_mat")
        cols_path = os.path.join(cols_dir, f"{cols_file_base}_cols.csv")
        identity_scores = read_csv(cols_path, skiprows=1).values.tolist()

        df = read_csv(
            matrix_path, delimiter=",", index_col=0, header=None, names=column_names
        )
        tick_text = df.index.tolist()
        update_document(doc_id, sequences_count=len(tick_text))
        data = df.to_numpy()

        diag_mask = eye(data.shape[0], dtype=bool)
        data_no_diag = where(diag_mask, nan, data)
        min_val = int(nanmin(data_no_diag))
        max_val = int(nanmax(data_no_diag))

        # TODO might be able to make one tick text object for both to use?
        return data, tick_text, min_val, max_val, identity_scores, stats_df

    def get_data(self, doc_id: str):
        data, tick_text, min_val, max_val, identity_scores, stats_df = (
            self.load_data_and_stats(doc_id)
        )
        heat_data = DataFrame(data, index=tick_text)
        parsedData = heat_data.values.tolist()

        data_to_dump = dict(
            metadata=dict(minVal=min_val, maxVal=max_val),
            data=([tick_text] + parsedData),
            identity_scores=identity_scores,
            full_stats=stats_df.values.tolist()
        )
        return json.dumps(data_to_dump)

    def new_doc(self):
        id = make_doc_id()
        new_document(id)
        return id

    def get_doc(self, doc_id: str):
        doc = get_document(doc_id)
        return doc._asdict() if doc else None

    def save_doc(self, doc_id: str, path: str, save_as: bool = False):
        print(doc_id, path, save_as)
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
                doc_id,
                filename=path,
                basename=basename,
                filetype="application/vnd.sdt"
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
            distribution=args["distribution"],
        )
        doc = get_document(args["id"])
        if doc == None:
            raise Exception(f"Expected to find document: {args['id']}")
        save_doc_settings(doc)

    def close_doc(self, doc_id: str):
        remove_document(doc_id)

    def set_window_title(self, title: str):
        assert_window().title = f"{title}{window_title_suffix}"

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
    #map(lambda doc: doc.cleanup(), get_state().documents)
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
            cores=multiprocessing.cpu_count(),
            memory=psutil.virtual_memory().total,
        ),
        on_update=lambda _: push_backend_state(window),
    )

    webview.start(debug=get_state().debug, private_mode=False)
