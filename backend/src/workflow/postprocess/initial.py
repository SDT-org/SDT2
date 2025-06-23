from pandas.core.frame import DataFrame
from export_utils import save_cols_to_csv, save_matrix_to_csv, save_stats_to_csv
from workflow.postprocess.utils import get_seq_stats

from workflow.models import RunSettings, WorkflowResult


def run(result: WorkflowResult, settings: RunSettings) -> WorkflowResult:
    doc_paths = settings.doc_paths
    matrix, ordered_ids = result.matrix, result.ordered_ids

    seq_stats = get_seq_stats(result.records)
    df = DataFrame(matrix, index=ordered_ids, columns=ordered_ids)
    save_cols_to_csv(df, doc_paths.columns)
    save_stats_to_csv(seq_stats, doc_paths.stats)
    save_matrix_to_csv(df, doc_paths.matrix, doc_paths.triangle)

    return result
