import subprocess
import os
from typing import Callable
from workflow.models import RunSettings, WorkflowResult

from transformations import lzani_tsv_to_distance_matrix


def run(
    result: WorkflowResult,
    settings: RunSettings,
    set_progress: Callable[[int], None],
    canceled,
) -> WorkflowResult:
    cmd = [
        os.path.abspath(os.path.normpath(settings.lzani.exec_path)),
        "all2all",
        "--in-fasta",
        settings.fasta_path,
        "--out",
        settings.doc_paths.lzani_results,
        "--out-ids",
        settings.doc_paths.lzani_results_ids,
        "--out-format",
        "complete",
        "--threads",
        "0",
        "--verbose",
        "2",
    ]
    #Not actual alignments, but can give insights into the analysis
    if settings.export_alignments and settings.alignment_export_path:
        alignment_file = os.path.join(settings.alignment_export_path, "lzani_alignments.txt")
        cmd.extend(["--out-alignment", alignment_file])

        os.makedirs(settings.alignment_export_path, exist_ok=True)
    
    # Add user entered parameters
    if settings.lzani.aw is not None:
        cmd.extend(["--aw", str(settings.lzani.aw)])
        
    if settings.lzani.am is not None:
        cmd.extend(["--am", str(settings.lzani.am)])
        
    if settings.lzani.mal is not None:
        cmd.extend(["--mal", str(settings.lzani.mal)])
        
    if settings.lzani.msl is not None:
        cmd.extend(["--msl", str(settings.lzani.msl)])
        
    if settings.lzani.mrd is not None:
        cmd.extend(["--mrd", str(settings.lzani.mrd)])
        
    if settings.lzani.mqd is not None:
        cmd.extend(["--mqd", str(settings.lzani.mqd)])
        
    if settings.lzani.reg is not None:
        cmd.extend(["--reg", str(settings.lzani.reg)])
        
    if settings.lzani.ar is not None:
        cmd.extend(["--ar", str(settings.lzani.ar)])
    
    process = subprocess.Popen(
        cmd,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
    )

    # Use communicate() to avoid pipe buffer issues
    try:
        stdout, stderr = process.communicate(timeout=300)  # 5 minute timeout
    except subprocess.TimeoutExpired:
        process.kill()
        stdout, stderr = process.communicate()
        return result._replace(error="LZ-ANI timed out after 5 minutes")

    if process.returncode != 0:
        return result._replace(error=f"Error running LZ-ANI: {stderr.strip()}")

    matrix, ids = lzani_tsv_to_distance_matrix(
        settings.doc_paths.lzani_results,
        settings.doc_paths.lzani_results_ids,
        settings.lzani.score_type,
    )

    return result._replace(distance_matrix=matrix, ordered_ids=ids)
