from multiprocessing.pool import Pool
from typing import Callable
import parasail
from functools import partial
from itertools import combinations_with_replacement as cwr
from numpy import zeros

from workflow.models import RunSettings, WorkflowResult
from workflow.analyze.scoring_matrices import SIMPLE_MATRICES


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

    matrix_to_use = settings.parasail.scoring_matrix or (
        "blosum62" if result.is_aa else "simple_2_-1"
    )
    open_penalty = settings.parasail.open_penalty or (10 if result.is_aa else 8)
    extend_penalty = settings.parasail.extend_penalty or 1

    # For alignment export, we need to use traceback regardless of SIMD support
    use_traceback = settings.export_alignments or not supports_striped_32()

    bound_process_pair = partial(
        process_pair,
        matrix_id=matrix_to_use,
        open_penalty=open_penalty,
        extend_penalty=extend_penalty,
        use_traceback=use_traceback,
        export_alignments=settings.export_alignments,
        alignment_export_path=settings.alignment_export_path,
    )

    # Create alignment export directory if needed
    if settings.export_alignments and settings.alignment_export_path:
        import os

        os.makedirs(settings.alignment_export_path, exist_ok=True)

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


def process_pair(
    id_sequence_pair,
    matrix_id,
    open_penalty,
    extend_penalty,
    use_traceback=False,
    export_alignments=False,
    alignment_export_path="",
):
    id1 = id_sequence_pair[0][0]
    id2 = id_sequence_pair[0][1]
    seq1 = id_sequence_pair[1][0]
    seq2 = id_sequence_pair[1][1]

    if id1 == id2:
        return id_sequence_pair[0], 0.0

    # Handle simple matrices
    if matrix_id in SIMPLE_MATRICES:
        matrix_info = SIMPLE_MATRICES[matrix_id]
        scoring_matrix = parasail.matrix_create(
            "ACGT", matrix_info["match"], matrix_info["mismatch"]
        )
    else:
        # Handle built-in parasail matrices
        matrix_map = {
            # BLOSUM matrices
            "blosum45": parasail.blosum45,
            "blosum62": parasail.blosum62,
            "blosum80": parasail.blosum80,
            "blosum90": parasail.blosum90,
            # PAM matrices
            "pam30": parasail.pam30,
            "pam70": parasail.pam70,
            "pam120": parasail.pam120,
            "pam250": parasail.pam250,
            # Nucleotide matrices
            "dnafull": parasail.dnafull,
            "nuc44": parasail.nuc44,
        }
        scoring_matrix = matrix_map.get(matrix_id, parasail.blosum62)

    # Use traceback if requested or if SIMD not supported
    if use_traceback:
        result = get_traceback_score(
            seq1, seq2, open_penalty, extend_penalty, scoring_matrix
        )

        # Export alignment if requested
        if export_alignments and alignment_export_path:
            try:
                alignment_result = parasail.nw_trace(
                    seq1, seq2, open_penalty, extend_penalty, scoring_matrix
                )
                import os

                filename = f"{id1}_vs_{id2}.txt".replace("/", "_").replace("\\", "_")
                filepath = os.path.join(alignment_export_path, filename)
                with open(filepath, "w") as f:
                    f.write(f">{id1}\n{alignment_result.traceback.query}\n")
                    f.write(f">{id2}\n{alignment_result.traceback.ref}\n")
            except Exception as e:
                print(f"Failed to export alignment for {id1} vs {id2}: {e}")
    else:
        result = get_stats_score(
            seq1, seq2, open_penalty, extend_penalty, scoring_matrix
        )

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
