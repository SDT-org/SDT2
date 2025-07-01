from typing import Dict
from Bio.SeqUtils import gc_fraction
from pandas.core.frame import DataFrame
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

    seq_stats = get_seq_stats(result.seq_dict)
    df = DataFrame(distance_matrix, index=ids_to_use, columns=ids_to_use)
    save_cols_to_csv(df, doc_paths.columns)
    save_stats_to_csv(seq_stats, doc_paths.stats)
    save_matrix_to_csv(df, doc_paths.matrix, doc_paths.triangle)
    save_seq_dict_to_json(result.seq_dict, doc_paths.seq_dict)

    return result


def get_seq_stats(seq_dict: Dict[str, str]) -> dict[str, list[float]]:
    seq_stats = {}
    for id, seq in seq_dict.items():
        gcCount = round(gc_fraction(seq), 2) * 100
        genLen = len(str(seq))
        seq_stats.setdefault(id, [])
        seq_stats[id].append(gcCount)
        seq_stats[id].append(genLen)
    return seq_stats
