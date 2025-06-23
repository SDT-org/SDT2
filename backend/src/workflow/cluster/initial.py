from workflow.models import RunSettings, WorkflowResult


def run(result: WorkflowResult, settings: RunSettings) -> WorkflowResult:
    cluster_method = settings.cluster_method
    final_ordered_ids, final_matrix_for_df = seq_ids_in_order, aln_scores
    ##3 lets bring this out as a compoennt for LZAZANI and parasail alike, just do some handling
    if cluster_method is not None and len(seq_ids_in_order) > 1:
        current_dist_matrix_for_clustering = np.copy(dist_scores)
        np.fill_diagonal(current_dist_matrix_for_clustering, 0)
        condensed_dist_matrix = [
            current_dist_matrix_for_clustering[i, j]
            for i in range(len(seq_ids_in_order))
            for j in range(i + 1, len(seq_ids_in_order))
        ]

        if condensed_dist_matrix:
            try:
                linked = hierarchy.linkage(
                    np.array(condensed_dist_matrix), method=cluster_method
                )
                dendro_data = hierarchy.dendrogram(
                    linked, orientation="right", no_plot=True
                )
                reordered_original_indices = [
                    int(i_str) for i_str in dendro_data["ivl"]
                ]
                final_ordered_ids = [
                    seq_ids_in_order[i] for i in reordered_original_indices
                ]
                final_matrix_for_df = aln_scores[
                    np.ix_(reordered_original_indices, reordered_original_indices)
                ]
            except Exception as e:
                print(
                    f"Warning: Clustering failed: {e}. Using original order.",
                    file=sys.stderr,
                )
                final_ordered_ids, final_matrix_for_df = seq_ids_in_order, aln_scores
        else:
            final_ordered_ids, final_matrix_for_df = seq_ids_in_order, aln_scores

    return result
