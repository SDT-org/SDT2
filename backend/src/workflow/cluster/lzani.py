from workflow.models import RunSettings, WorkflowResult
from scipy.cluster import hierarchy


def run(result: WorkflowResult, settings: RunSettings) -> WorkflowResult:
    matrix = result.matrix
    dist_matrix = 100.0 - matrix.values
    np.fill_diagonal(dist_matrix, 0)
    condensed_dist = [
        dist_matrix[i, j] for i in range(len(matrix)) for j in range(i + 1, len(matrix))
    ]
    ##--need to plug in the reorder logic from parasail here as a component
    linked = hierarchy.linkage(np.array(condensed_dist), method=settings.cluster_method)
    dendro_data = hierarchy.dendrogram(linked, orientation="right", no_plot=True)
    reordered_indices = [int(i) for i in dendro_data["ivl"]]
    seq_ids = matrix.index.tolist()
    reordered_ids = [seq_ids[i] for i in reordered_indices]
    matrix = matrix.loc[reordered_ids, reordered_ids]

    return result
