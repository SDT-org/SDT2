from datetime import datetime
from multiprocessing.sharedctypes import Synchronized
import os
import platform
import time
import numpy
from pandas import DataFrame
import psutil

from config import app_version
from export_utils import save_run_settings_to_json
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
    return parse.run(result, fasta_path)


def run_process(
    workflow_run: WorkflowRun, cancel_event: Synchronized
) -> WorkflowResult:
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

    # Skip clustering and post-processing for large datasets in lzani and vclust workflows
    sequence_count = len(result.ordered_ids)
    if settings.analysis_method in ["lzani", "vclust"] and sequence_count > 1000:
        print(f"Skipping clustering and post-processing for large dataset ({sequence_count} sequences > 1000)")
        
        # For large datasets, the only valid view is UMAP.
        # We must update the document state to switch the view automatically.
        from app_state import update_document
        # Extract doc_id from the document path (it's the parent directory name)
        doc_id = os.path.basename(os.path.dirname(settings.doc_paths.matrix))
        update_document(doc_id, dataView='umap')
        
        # We must also skip the finalization step for large datasets
        workflow_run.set_stage("Finalizing")
        workflow_run.set_progress(100)
        
        end_time, end_counter = datetime.now(), time.perf_counter()

        # Build comprehensive LZ-ANI settings
        lzani_settings = {
            "score_type": settings.lzani.score_type,
        }
        if settings.lzani.aw is not None:
            lzani_settings["aw"] = settings.lzani.aw
        if settings.lzani.am is not None:
            lzani_settings["am"] = settings.lzani.am
        if settings.lzani.mal is not None:
            lzani_settings["mal"] = settings.lzani.mal
        if settings.lzani.msl is not None:
            lzani_settings["msl"] = settings.lzani.msl
        if settings.lzani.mrd is not None:
            lzani_settings["mrd"] = str(settings.lzani.mrd)
        if settings.lzani.mqd is not None:
            lzani_settings["mqd"] = str(settings.lzani.mqd)
        if settings.lzani.reg is not None:
            lzani_settings["reg"] = settings.lzani.reg
        if settings.lzani.ar is not None:
            lzani_settings["ar"] = settings.lzani.ar

        run_settings = {
            "analysis_method": settings.analysis_method,
            "cluster_method": settings.cluster_method,
            "lzani": lzani_settings,
        }
        
        # Add vclust settings if applicable
        if settings.analysis_method == "vclust" and settings.vclust:
            run_settings["vclust"] = {
                "kmer_min_similarity": settings.vclust.kmer_min_similarity,
                "kmer_min_kmers": settings.vclust.kmer_min_kmers,
                "kmer_fraction": settings.vclust.kmer_fraction,
                "cdhit_threshold": settings.vclust.cdhit_threshold,
            }

        save_run_settings_to_json(run_settings, settings.doc_paths.run_settings)

        file_name = os.path.basename(settings.fasta_path)
        summary_text = output_summary(
            file_name, start_time, end_time, start_counter, end_counter, settings, result
        )
        with open(settings.doc_paths.summary, "w") as file:
            file.write(summary_text)
        print(f"Elapsed time: {friendly_total_time(end_counter - start_counter)}")
        
        # For vclust, store the sparse matrix in cache to avoid JSON serialization
        if settings.analysis_method == "vclust" and result.distance_matrix is not None:
            from workflow.vclust_cache import store_vclust_matrix
            store_vclust_matrix(doc_id, result.distance_matrix, result.ordered_ids)
        
        return result

    elif settings.cluster_method and settings.cluster_method != "None":
        workflow_run.set_stage("Clustering")
        workflow_run.set_progress(None)

        result = cluster.run(result, settings, workflow_run.set_progress)
        if result.error:
            return result

    workflow_run.set_stage("Finalizing")
    workflow_run.set_progress(100)

    result = postprocess.run(result, settings)
    if result.error:
        return result

    end_time, end_counter = datetime.now(), time.perf_counter()

    # Build comprehensive LZ-ANI settings
    lzani_settings = {
        "score_type": settings.lzani.score_type,
    }
    if settings.lzani.aw is not None:
        lzani_settings["aw"] = settings.lzani.aw
    if settings.lzani.am is not None:
        lzani_settings["am"] = settings.lzani.am
    if settings.lzani.mal is not None:
        lzani_settings["mal"] = settings.lzani.mal
    if settings.lzani.msl is not None:
        lzani_settings["msl"] = settings.lzani.msl
    if settings.lzani.mrd is not None:
        lzani_settings["mrd"] = str(settings.lzani.mrd)
    if settings.lzani.mqd is not None:
        lzani_settings["mqd"] = str(settings.lzani.mqd)
    if settings.lzani.reg is not None:
        lzani_settings["reg"] = settings.lzani.reg
    if settings.lzani.ar is not None:
        lzani_settings["ar"] = settings.lzani.ar

    run_settings = {
        "analysis_method": settings.analysis_method,
        "cluster_method": settings.cluster_method,
        "lzani": lzani_settings,
    }
    
    # Add vclust settings if applicable
    if settings.analysis_method == "vclust" and settings.vclust:
        run_settings["vclust"] = {
            "kmer_min_similarity": settings.vclust.kmer_min_similarity,
            "kmer_min_kmers": settings.vclust.kmer_min_kmers,
            "kmer_fraction": settings.vclust.kmer_fraction,
            "cdhit_threshold": settings.vclust.cdhit_threshold,
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


def output_summary(
    file_name, start_time, end_time, start_counter, end_counter, settings, result
):
    build_type = f"{platform.system()} {platform.release()} {platform.machine()}"
    total_cores, total_ram = os.cpu_count(), psutil.virtual_memory().total / (1024**3)
    total_counter = end_counter - start_counter

    if settings.analysis_method == "parasail":

        if settings.parasail.scoring_matrix:
            scoring_matrix = settings.parasail.scoring_matrix
        else:
            scoring_matrix = "BLOSUM62" if result.is_aa else "1-1"

        open_penalty = settings.parasail.open_penalty or 10
        extend_penalty = settings.parasail.extend_penalty or 1
        aligner_name = "PARASAIL"
        aligner_params = f"Needleman-Wunsch (stats) algorithm, Scoring matrix: {scoring_matrix}, Open penalty: {open_penalty}, Extend penalty: {extend_penalty}"
    elif settings.analysis_method == "lzani":
        score_type = settings.lzani.score_type.upper()

        aw = settings.lzani.aw if settings.lzani.aw is not None else 5
        am = settings.lzani.am if settings.lzani.am is not None else 2
        mal = settings.lzani.mal if settings.lzani.mal is not None else 50
        msl = settings.lzani.msl if settings.lzani.msl is not None else 50
        mrd = settings.lzani.mrd if settings.lzani.mrd is not None else 0.1
        mqd = settings.lzani.mqd if settings.lzani.mqd is not None else 0.01
        reg = settings.lzani.reg if settings.lzani.reg is not None else 0
        ar = settings.lzani.ar if settings.lzani.ar is not None else 0.95

        aligner_name = "LZ-ANI"
        aligner_params = f"Score type: {score_type}, aw={aw}, am={am}, mal={mal}, msl={msl}, mrd={mrd}, mqd={mqd}, reg={reg}, ar={ar}"
    elif settings.analysis_method == "vclust":
        aligner_name = "VCLUST"
        vclust_params = []
        if settings.vclust:
            vclust_params.append(f"K-mer min similarity: {settings.vclust.kmer_min_similarity}")
            vclust_params.append(f"K-mer min kmers: {settings.vclust.kmer_min_kmers}")
            vclust_params.append(f"K-mer fraction: {settings.vclust.kmer_fraction}")
            vclust_params.append(f"CD-HIT threshold: {settings.vclust.cdhit_threshold}")
        
        # Also include LZ-ANI parameters used in vclust
        score_type = settings.lzani.score_type.upper()
        aw = settings.lzani.aw if settings.lzani.aw is not None else 5
        am = settings.lzani.am if settings.lzani.am is not None else 2
        mal = settings.lzani.mal if settings.lzani.mal is not None else 50
        msl = settings.lzani.msl if settings.lzani.msl is not None else 50
        mrd = settings.lzani.mrd if settings.lzani.mrd is not None else 0.1
        mqd = settings.lzani.mqd if settings.lzani.mqd is not None else 0.01
        reg = settings.lzani.reg if settings.lzani.reg is not None else 0
        ar = settings.lzani.ar if settings.lzani.ar is not None else 0.95
        
        vclust_params.append(f"LZ-ANI score type: {score_type}, aw={aw}, am={am}, mal={mal}, msl={msl}, mrd={mrd}, mqd={mqd}, reg={reg}, ar={ar}")
        aligner_params = ", ".join(vclust_params)
    else:
        aligner_name = "UNKNOWN"
        aligner_params = f"Unknown analysis method: {settings.analysis_method}"

    return f"""
    SDT {app_version} release for {platform.system()}
    Developed by Michael Lund and Josiah Ivey; Brejnev Muhire,
    Darren Martin, Simona Kraberger, Qiyun Zhu, Pierre Lefeuvre, Jean-Michele Lett, Philippe Roumagnac, Arvind Varsani
    System info: Host: {platform.node()}, OS: {build_type}, CPU: {platform.processor()} - {total_cores} cores, Memory: {total_ram:.2f} GB RAM
    Run info for {file_name}: Start: {start_time.strftime("%b %d %Y, %I:%M %p %Z")}, End: {end_time.strftime("%b %d %Y, %I:%M %p %Z")}, Total: {friendly_total_time(total_counter)}
    Using the {aligner_name} aligner. Run parameters: {aligner_params}
    """
