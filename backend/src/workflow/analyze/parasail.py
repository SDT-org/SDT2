from multiprocessing.pool import Pool
from typing import Callable
import parasail
import sys
from functools import partial
from itertools import combinations_with_replacement as cwr
from numpy import zeros

from workflow.models import RunSettings, WorkflowResult


def run(
    result: WorkflowResult,
    settings: RunSettings,
    set_progress: Callable[[int], None],
    canceled,
) -> WorkflowResult:
    seq_ids = list(result.seq_dict.keys())
    seq_dict = result.seq_dict
    n = len(seq_ids)
    dist_scores = zeros((n, n))
    order = {seq_id: i for i, seq_id in enumerate(seq_ids)}
    combos = list(cwr(seq_ids, 2))
    id_sequence_pairs = []
    min_score = 0

    for ids_tuple in combos:
        id_sequence_pairs.append(
            [ids_tuple, [seq_dict[ids_tuple[0]], seq_dict[ids_tuple[1]]]]
        )

    total_pairs = len(id_sequence_pairs)

    print(f"\rNumber of sequences: {len(seq_ids)}\r", flush=True)
    print(f"\rNumber of pairs: {total_pairs}\r", flush=True)
    bound_process_pair = partial(process_pair, is_aa=result.is_aa)

    with Pool(settings.parasail.process_count) as pool:
        results = pool.imap(bound_process_pair, id_sequence_pairs)
        for i, pair_result in enumerate(results, 1):
            if canceled.value:
                pool.close()
                pool.terminate()
                pool.join()
                canceled.value = False
                return result._replace(distance_matrix=None, error="PROCESS_CANCELED")
            [seqid1, seqid2], score = pair_result
            seqid1, seqid2 = order[seqid1], order[seqid2]
            dist_scores[seqid1, seqid2] = score
            dist_scores[seqid2, seqid1] = score
            if min_score is None or score < min_score:
                min_score = score
            progress = int((i / total_pairs) * 100)
            set_progress(progress)

    return result._replace(
        ordered_ids=list(order.keys()), distance_matrix=dist_scores, min_score=min_score
    )


def supports_striped_32():
    return (
        parasail.can_use_sse41() or parasail.can_use_avx2() or parasail.can_use_neon()
    )


def process_pair(id_sequence_pair, is_aa):
    id1 = id_sequence_pair[0][0]
    id2 = id_sequence_pair[0][1]
    seq1 = id_sequence_pair[1][0]
    seq2 = id_sequence_pair[1][1]

    if id1 == id2:
        return id_sequence_pair[0], 0.0

    open_penalty = 10 if is_aa else 8
    extend_penalty = 1
    matrix_to_use = parasail.blosum62 if is_aa else parasail.pam250
    alignment_function = (
        get_stats_score if supports_striped_32() else get_traceback_score
    )
    result = alignment_function(seq1, seq2, open_penalty, extend_penalty, matrix_to_use)

    return id_sequence_pair[0], result


def get_stats_score(seq1, seq2, open_penalty, extend_penalty, matrix) -> float:
    # Use 32-bit version to prevent score overflow with long sequences
    result = parasail.nw_stats_striped_32(
        seq1, seq2, open_penalty, extend_penalty, matrix
    )
    num_ungapped_cols = len(seq1) + len(seq2) - result.length
    identity_score_percent = (float(result.matches) / num_ungapped_cols) * 100.0
    distance_score = 100.0 - identity_score_percent
    return distance_score


def get_similarity(seq1, seq2):
    if len(seq1) != len(seq2):
        raise ValueError("Strings must be of equal length.")
    dist = 0
    gaps = 0
    alns = zip(seq1, seq2)
    for a, b in alns:
        if a != "-" and b != "-":
            if a != b:
                dist += 1
        else:
            gaps += 1

    similarity = float((float(dist)) / (len(seq1) - gaps))

    return similarity


def get_traceback_score(seq1, seq2, open_penalty, extend_penalty, matrix) -> float:
    try:
        result = parasail.nw_trace(seq1, seq2, open_penalty, extend_penalty, matrix)
        query = result.traceback.query
    except:
        raise Exception("PARASAIL_TRACEBACK")

    similarity = get_similarity(query, result.traceback.ref)

    similarity_percent = similarity * 100.0
    return similarity_percent
