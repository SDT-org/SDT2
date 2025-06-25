import subprocess
import os
import time
from typing import Callable
import pandas as pd
import numpy as np
from workflow.models import RunSettings, WorkflowResult
from pandas import DataFrame


def run(
    result: WorkflowResult,
    settings: RunSettings,
    set_progress: Callable[[int], None],
    pool,
    cancel_event,
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

    while process.poll() is None:
        if cancel_event.is_set():
            process.terminate()
            try:
                process.wait(timeout=5)
            except subprocess.TimeoutExpired:
                process.kill()
            return result
        time.sleep(0.1)

    _, stderr = process.communicate()

    if process.returncode != 0:
        result.errors.append(f"Error running LZ-ANI: {stderr.strip()}")
        return result

    matrix, ids = lzani_to_full_matrix(
        settings.doc_paths.lzani_results,
        settings.doc_paths.lzani_results_ids,
        settings.lzani.score_type,
    )
    
    return result._replace(distance_matrix=matrix, ordered_ids=ids)


def lzani_to_full_matrix(results_tsv_path, ids_tsv_path, score_column="ani"):
    all_ids = pd.read_csv(ids_tsv_path, sep="\t")["id"].tolist()
    df = pd.read_csv(results_tsv_path, sep="\t")
    
    matrix = (
        df.pivot_table(
            index="query", columns="reference", values=score_column, aggfunc="first"
        )
        * 100
    )
    matrix = matrix.reindex(index=all_ids, columns=all_ids)
    matrix = matrix.fillna(25.0)
    matrix = matrix.replace(0.0, 25.0)
    
    symmetric_matrix = (matrix + matrix.T) / 2
    np.fill_diagonal(symmetric_matrix.values, 100.0)
    #returning np arrau and ids as tuple instead of DF
    return symmetric_matrix.values, all_ids
