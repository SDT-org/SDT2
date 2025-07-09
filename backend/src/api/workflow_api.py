import os
import sys
import platform
import psutil
from multiprocessing import Manager, cpu_count as get_cpu_count

sys.path.append(os.path.join(os.path.dirname(__file__), ".."))

from workflow.models import (
    LzaniSettings,
    ParasailSettings,
    RunSettings,
    WorkflowRun,
    WorkflowResult,
)
from workflow.runner import run_process
from config.paths import build_document_paths
from config.constants import default_window_title

is_windows = platform.system() == "Windows"
try:
    cpu_count = get_cpu_count()
except:
    cpu_count = 1


def get_lzani_exec_path():
    current_script_dir = os.path.dirname(os.path.abspath(__file__))
    backend_src_dir = os.path.dirname(current_script_dir)
    backend_dir = os.path.dirname(backend_src_dir)
    bin_dir = os.path.join(backend_dir, "bin")
    lzani_executable_name = "lz-ani.exe" if is_windows else "lz-ani"
    lzani_executable_path = os.path.join(bin_dir, lzani_executable_name)
    if not os.path.exists(lzani_executable_path):
        raise FileNotFoundError(
            f"LZ-ANI executable not found at {lzani_executable_path}"
        )
    return lzani_executable_path


def get_compute_stats(workflow_result: WorkflowResult, get_state):
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


def do_cancel_run(cancel_event, get_state, update_document, set_state, workflow_runs):
    if cancel_event:
        try:
            cancel_event.value = True
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


def start_workflow_run(args: dict, get_state, get_document, update_document,
                      set_state, parsed_workflow_results, workflow_runs):
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
                return
            raise Exception(f"Workflow processing step failed: {result.error}")

    set_state(active_run_document_id=None)
    del workflow_runs[doc.id]
    update_document(doc.id, view="viewer")


def get_workflow_run_status(doc_id: str, workflow_runs):
    workflow_run = workflow_runs.get(doc_id)
    if not workflow_run:
        raise Exception(f"No workflow run found for document ID: {doc_id}")
    return {"stage": workflow_run.stage, "progress": workflow_run.progress}


def cancel_run(doc_id: str, run_settings: str, cancel_event, get_state, 
               update_document, set_state, workflow_runs, reset_state, window):
    if cancel_event:
        try:
            cancel_event.value = True
        except:
            pass
    doc_id_active = get_state().active_run_document_id
    if doc_id_active:
        del workflow_runs[doc_id_active]
        update_document(
            doc_id_active,
            view="runner",
            doc_id=doc_id_active,
            progress=0,
            pair_progress=0,
            pair_count=0,
        )
    set_state(active_run_document_id=None)
    print("Run canceled")
    
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
        if window:
            window.title = default_window_title
