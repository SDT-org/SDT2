import parasail


def get_stats_score(seq1, seq2, open_penalty, extend_penalty, matrix) -> float:
    # Use 32-bit version to prevent score overflow with long sequences
    result = parasail.nw_stats_striped_32(seq1, seq2, open_penalty, extend_penalty, matrix)
    num_ungapped_cols = len(seq1) + len(seq2) - result.length
    identity_score_percent = (float(result.matches) / num_ungapped_cols) * 100.0
    distance_score = 100.0 - identity_score_percent

    return distance_score