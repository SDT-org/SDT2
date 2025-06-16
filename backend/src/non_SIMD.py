import parasail

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
    # convert to percentile
    similarity_percentile = similarity * 100
    return similarity_percentile

def get_traceback_score(seq1, seq2, open_penalty, extend_penalty, matrix) -> float:
    try:
        result = parasail.nw_trace(seq1, seq2, open_penalty, extend_penalty, matrix)
        query = result.traceback.query
    except:
        raise Exception("PARASAIL_TRACEBACK")

    return get_similarity(query, result.traceback.ref)