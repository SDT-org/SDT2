import subprocess
import os
import time
from typing import Callable
import pandas as pd
import numpy as np
from workflow.models import RunSettings, WorkflowResult
from pandas import DataFrame
from transformations import lzani_tsv_to_distance_matrix


def run(
    result: WorkflowResult,
    settings: RunSettings,
    set_progress: Callable[[int], None],
    canceled,
) -> WorkflowResult:
    process = subprocess.Popen(
        [
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
            "0",  # 0 will auto-detect
            "--verbose",
            "2",
        ],
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
