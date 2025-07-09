from datetime import datetime
from multiprocessing.pool import Pool
import os
import platform
import time
import numpy
from pandas import DataFrame
import psutil

from config import app_version
from file_io.export_manager import save_run_settings_to_json
from utils import friendly_total_time
from workflow import analyze, cluster, parse, postprocess
from workflow.models import WorkflowResult, WorkflowRun


def run_parse(fasta_path: str) -> WorkflowResult:
    result = WorkflowResult(
        seq_dict={},
        ordered_ids=[],
        reordered_ids=[],
        max_sequence_length=0,
        warnings=[],
        error=None,
        distance_matrix=numpy.ndarray([]),
        similarity_matrix=DataFrame(),
        is_aa=None,
        min_score=0,
    )
    result = parse.run(result, fasta_path)
    if result.error:
        return result

    return result


def run_process(workflow_run: WorkflowRun, cancel_event) -> WorkflowResult:
    result = workflow_run.result
    settings = workflow_run.settings

    workflow_run.set_stage("Analyzing")
    workflow_run.set_progress(0)
    start_time, start_counter = (
        workflow_run.analyze_start_time or datetime.now(),
        workflow_run.analyze_start_counter or time.perf_counter(),
    )

    result = analyze.jobs[settings.analysis_method].run(
        result,
        settings,
        workflow_run.set_progress,
        cancel_event,
    )

    if result.error:
        return result

    if settings.cluster_method and settings.cluster_method != "None":
        workflow_run.set_stage("Clustering")
        workflow_run.set_progress(None)

        result = cluster.run(result, settings)
        if result.error:
            return result

    workflow_run.set_stage("Finalizing")
    workflow_run.set_progress(100)

    result = postprocess.run(result, settings)
    if result.error:
        return result

    end_time, end_counter = datetime.now(), time.perf_counter()

    run_settings = {
        "analysis_method": settings.analysis_method,
        "cluster_method": settings.cluster_method,
        "lzani": {
            "score_type": settings.lzani.score_type,
        },
    }

    save_run_settings_to_json(run_settings, settings.doc_paths.run_settings)

    file_name = os.path.basename(settings.fasta_path)
    summary_text = output_summary(
        file_name, start_time, end_time, start_counter, end_counter
    )
    with open(settings.doc_paths.summary, "w") as file:
        file.write(summary_text)
    print(f"Elapsed time: {friendly_total_time(end_counter - start_counter)}")

    return result


def output_summary(file_name, start_time, end_time, start_counter, end_counter):
    build_type = f"{platform.system()} {platform.release()} {platform.machine()}"
    total_cores, total_ram = os.cpu_count(), psutil.virtual_memory().total / (1024**3)
    total_counter = end_counter - start_counter
    return f"""
    SDT {app_version} release for {platform.system()}
    Developed by Michael Lund and Josiah Ivey; Brejnev Muhire,
    Darren Martin, Simona Kraberger, Qiyun Zhu, Pierre Lefeuvre, Jean-Michele Lett, Philippe Roumagnac, Arvind Varsani
    System info: Host: {platform.node()}, OS: {build_type}, CPU: {platform.processor()} - {total_cores} cores, Memory: {total_ram:.2f} GB RAM
    Run info for {file_name}: Start: {start_time.strftime("%b %d %Y, %I:%M %p %Z")}, End: {end_time.strftime("%b %d %Y, %I:%M %p %Z")}, Total: {friendly_total_time(total_counter)}
    Parasail: Using Needleman-Wunsch (stats) algorithm. Nucleotide: Open={13}, Extend={1} (BLOSUM62). Amino acid: Open={10}, Extend={1} (BLOSUM62).
    """
