from datetime import datetime
from multiprocessing.pool import Pool
import os
import platform
import time
import numpy

import psutil

from config import app_version
from workflow import analyze, cluster, parse, postprocess
from workflow.models import WorkflowResult, WorkflowRun


def run_parse(fasta_path: str) -> WorkflowResult:
    result = WorkflowResult(
        seq_dict={},
        ordered_ids=[],
        max_sequence_length=0,
        warnings=[],
        errors=[],
        matrix=numpy.ndarray([]),
        is_aa=False,
    )
    result = parse.run(result, fasta_path)
    if result.errors:
        return result
    print(result.seq_dict)
    return result


def run_process(workflow_run: WorkflowRun, pool: Pool, cancel_event) -> WorkflowResult:
    result = workflow_run.result
    settings = workflow_run.settings

    workflow_run.set_stage("Analyzing")
    start_time, start_counter = (
        workflow_run.analyze_start_time or datetime.now(),
        workflow_run.analyze_start_counter or time.perf_counter(),
    )

    result = analyze.jobs[settings.analysis_method].run(
        result,
        settings,
        workflow_run.set_progress,
        pool,
        cancel_event,
    )
    if result.errors:
        return result

    # Skip clustering for now...
    # if settings.cluster_method and settings.cluster_method != "None":
    #     workflow_run.set_stage("Clustering")

    #     result = cluster.initial.run(result, settings)
    #     if result.errors:
    #         return result

    workflow_run.set_stage("Finalizing")

    result = postprocess.run(result, settings)
    if result.errors:
        return result

    end_time, end_counter = datetime.now(), time.perf_counter()
    file_name = os.path.basename(settings.fasta_path)
    summary_text = output_summary(
        file_name, start_time, end_time, start_counter, end_counter
    )
    with open(settings.doc_paths.summary, "w") as file:
        file.write(summary_text)
    print(f"Elapsed time: {friendly_total_time(end_counter - start_counter)}")

    return result


def friendly_total_time(total_time):
    m, s = divmod(total_time, 60)
    return f"{int(m)}m {s:.2f}s" if m > 0 else f"{s:.2f}s"


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


# def fasta_alignments(seq_records, fname):
#     SeqIO.write(seq_records, fname, "fasta")
