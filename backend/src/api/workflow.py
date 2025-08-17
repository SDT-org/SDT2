from multiprocessing import Manager, cpu_count as get_cpu_count
import os
import sys

import psutil
from app import do_cancel_run
from document_paths import build_document_paths
from workflow.models import (
    LzaniSettings,
    ParasailSettings,
    VclustSettings,
    RunSettings,
    WorkflowResult,
    WorkflowRun,
)
from workflow.runner import run_process
from app_state import get_state, set_state, get_document, update_document
from app_globals import (
    get_workflow_runs,
    get_parsed_workflow_results,
    set_canceled,
)
from config import is_windows, is_compiled


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


def get_lzani_exec_path():  # TODO: move to dynamic config? or may need to be configurable by user
    lzani_executable_name = "lz-ani.exe" if is_windows else "lz-ani"

    if is_compiled:
        # When compiled, the executable is in the same directory as the main exe
        # but we need to go to backend/bin relative to the exe location
        exe_dir = os.path.dirname(sys.executable)
        lzani_executable_path = os.path.join(
            exe_dir, "backend", "bin", lzani_executable_name
        )
    else:
        # When running from source
        current_script_dir = os.path.dirname(
            os.path.abspath(__file__)
        )  # .../backend/src
        workflow_dir = os.path.dirname(current_script_dir)  # .../backend
        bin_dir = os.path.join(workflow_dir, "..", "bin")  # .../backend/bin
        lzani_executable_path = os.path.join(bin_dir, lzani_executable_name)

    if not os.path.exists(lzani_executable_path):
        raise FileNotFoundError(
            f"LZ-ANI executable not found at {lzani_executable_path}"
        )
    return lzani_executable_path


class Workflow:
    def start_workflow_run(self, args: dict):
        app_state = get_state()
        if app_state.debug:
            print("\nAPI args:", args)

        if app_state.active_run_document_id:
            raise Exception("Multiple runs are not supported")

        doc = get_document(args["doc_id"])

        parsed_result = get_parsed_workflow_results().get(doc.id)
        if not parsed_result:
            raise Exception(
                "Parsed workflow result not found. Close the document and reselect the file."
            )

        if parsed_result.error:
            raise Exception("Workflow has an error that must be resolved.")

        # Create vclust settings if it's a vclust run
        vclust_settings = None
        if args.get("analysisMethod") == "vclust":
            vclust_settings = VclustSettings(
                kmer_min_similarity=args.get("kmer_min_similarity", 0.30),
                kmer_min_kmers=args.get("kmer_min_kmers", 2),
                kmer_fraction=args.get("kmer_fraction", 0.5),
                cdhit_threshold=args.get("cdhit_threshold", 0.70),
            )

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
            vclust=vclust_settings,
            export_alignments=args.get("export_alignments", False),
            alignment_export_path=args.get("alignment_export_path", ""),
        )

        workflow_run = WorkflowRun(
            result=parsed_result,
            settings=settings,
            progress=0,
        )

        get_workflow_runs()[doc.id] = workflow_run

        update_document(doc.id, view="loader")
        with Manager() as manager:
            canceled = manager.Value("b", False)
            set_canceled(canceled)
            result = run_process(workflow_run, canceled)

            if result.error:
                if result.error == "PROCESS_CANCELED":
                    return
                raise Exception(f"Workflow processing step failed: {result.error}")

        set_state(active_run_document_id=None)
        workflow_runs = get_workflow_runs()
        if doc.id in workflow_runs:
            del workflow_runs[doc.id]
        update_document(doc.id, view="viewer")

    def get_workflow_run_status(self, doc_id: str):
        workflow_run = get_workflow_runs().get(doc_id)
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
            from api.documents import Documents

            Documents().reset_state()
