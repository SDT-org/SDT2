import pandas as pd
from typing import Dict
from Bio.SeqUtils import gc_fraction
from pandas.core.frame import DataFrame
from scipy.sparse import issparse
from export_utils import (
    save_cols_to_csv,
    save_matrix_to_csv,
    save_stats_to_csv,
    save_seq_dict_to_json,
)
from workflow.models import RunSettings, WorkflowResult


def run(result: WorkflowResult, settings: RunSettings) -> WorkflowResult:
    doc_paths = settings.doc_paths

    # Use reordered_ids if clustering was performed, otherwise use ordered_ids
    ids_to_use = result.reordered_ids if result.reordered_ids else result.ordered_ids
    distance_matrix = result.distance_matrix

    seq_stats = get_seq_stats(result.seq_dict, result.is_aa or False)
    
    # Check if this is a large dataset with dummy matrix (1x1)
    sequence_count = len(ids_to_use)
    if sequence_count > 2500 and distance_matrix.shape == (1, 1):
        print(f"Skipping matrix save for large dataset ({sequence_count} sequences > 2500)")
        # Create minimal files for compatibility
        # Save empty columns file
        with open(doc_paths.columns, 'w') as f:
            f.write("")
        # Save stats as usual
        save_stats_to_csv(seq_stats, doc_paths.stats)
        # Save dummy matrix files to avoid errors
        with open(doc_paths.matrix, 'w') as f:
            f.write("")
        with open(doc_paths.triangle, 'w') as f:
            f.write("")
    else:
        # For the vclust workflow (which produces a sparse matrix), we skip saving the matrix to disk.
        if issparse(distance_matrix):
            print("Sparse matrix detected (vclust workflow), skipping matrix save.")
            save_stats_to_csv(seq_stats, doc_paths.stats)
            with open(doc_paths.columns, 'w') as f:
                f.write("")
        # For other workflows, we process the dense matrix as usual.
        else:
            df = DataFrame(distance_matrix, index=ids_to_use, columns=ids_to_use)
            save_cols_to_csv(df, doc_paths.columns)
            save_stats_to_csv(seq_stats, doc_paths.stats)
            save_matrix_to_csv(df, doc_paths.matrix, doc_paths.triangle)
    
    save_seq_dict_to_json(result.seq_dict, doc_paths.seq_dict)

    return result


def get_seq_stats(seq_dict: Dict[str, str], is_aa: bool) -> dict[str, list[float]]:
    seq_stats = {}
    for id, seq in seq_dict.items():
        if is_aa:
            gcCount = 0.0
        else:
            gcCount = round(gc_fraction(seq), 2) * 100
        genLen = len(str(seq))
        seq_stats.setdefault(id, [])
        seq_stats[id].append(gcCount)
        seq_stats[id].append(genLen)
    return seq_stats
