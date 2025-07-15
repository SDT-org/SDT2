from datetime import datetime
from multiprocessing.pool import Pool
import os
import platform
import time
import numpy
from pandas import DataFrame
import psutil

from config.config import app_version
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
        file_name, start_time, end_time, start_counter, end_counter, settings, result
    )
    with open(settings.doc_paths.summary, "w") as file:
        file.write(summary_text)
    print(f"Elapsed time: {friendly_total_time(end_counter - start_counter)}")

    return result


def output_summary(file_name, start_time, end_time, start_counter, end_counter, settings, result):
    build_type = f"{platform.system()} {platform.release()} {platform.machine()}"
    total_cores, total_ram = os.cpu_count(), psutil.virtual_memory().total / (1024**3)
    total_counter = end_counter - start_counter
    
    # Build analysis method details
    analysis_details = ""
    if settings.analysis_method == "parasail":
        # Determine if using custom or default settings
        if settings.parasail.scoring_matrix is not None:
            matrix = settings.parasail.scoring_matrix
        else:
            matrix = "BLOSUM62"  # default
        
        if settings.parasail.open_penalty is not None:
            open_penalty = settings.parasail.open_penalty
        else:
            open_penalty = 10 if result.is_aa else 13  # defaults
            
        if settings.parasail.extend_penalty is not None:
            extend_penalty = settings.parasail.extend_penalty
        else:
            extend_penalty = 1  # default
            
        analysis_details = f"""
    Parasail: Using Needleman-Wunsch (stats) algorithm
    - Scoring Matrix: {matrix}
    - Open Penalty: {open_penalty}
    - Extend Penalty: {extend_penalty}
    - Process Count: {settings.parasail.process_count}
    - Sequence Type: {'Amino Acid' if result.is_aa else 'Nucleotide'}"""
    
    elif settings.analysis_method == "lzani":
        # Build LZANI parameters string
        params = []
        if settings.lzani.aw is not None:
            params.append(f"aw={settings.lzani.aw}")
        if settings.lzani.am is not None:
            params.append(f"am={settings.lzani.am}")
        if settings.lzani.mal is not None:
            params.append(f"mal={settings.lzani.mal}")
        if settings.lzani.msl is not None:
            params.append(f"msl={settings.lzani.msl}")
        if settings.lzani.mrd is not None:
            params.append(f"mrd={settings.lzani.mrd}")
        if settings.lzani.mqd is not None:
            params.append(f"mqd={settings.lzani.mqd}")
        if settings.lzani.reg is not None:
            params.append(f"reg={settings.lzani.reg}")
        if settings.lzani.ar is not None:
            params.append(f"ar={settings.lzani.ar}")
            
        params_str = "default" if not params else ", ".join(params)
        
        analysis_details = f"""
    LZANI: Using {settings.lzani.score_type.upper()} score type
    - Parameters: {params_str}"""
    
    # Clustering details
    cluster_details = ""
    if settings.cluster_method and settings.cluster_method != "None":
        cluster_details = f"\n    Clustering: {settings.cluster_method}"
    
    return f"""
    SDT {app_version} release for {platform.system()}
    Developed by Michael Lund and Josiah Ivey; Brejnev Muhire,
    Darren Martin, Simona Kraberger, Qiyun Zhu, Pierre Lefeuvre, Jean-Michele Lett, Philippe Roumagnac, Arvind Varsani
    System info: Host: {platform.node()}, OS: {build_type}, CPU: {platform.processor()} - {total_cores} cores, Memory: {total_ram:.2f} GB RAM
    Run info for {file_name}: Start: {start_time.strftime("%b %d %Y, %I:%M %p %Z")}, End: {end_time.strftime("%b %d %Y, %I:%M %p %Z")}, Total: {friendly_total_time(total_counter)}{analysis_details}{cluster_details}
    """
