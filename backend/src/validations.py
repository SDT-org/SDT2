from Bio import SeqIO


def validate_fasta(filepath, filetype=None):
    if filetype and filetype != "text/fasta":
        return False, "INVALID_FILE_TYPE"
    try:
        records = list(SeqIO.parse(filepath, "fasta"))
    except Exception:
        return False, "INVALID_FASTA_FORMAT"

    if len(records) < 2:
        return False, "NOT_ENOUGH_SEQUENCES"

    sequence_names = set()
    for record in records:
        if not record.seq:
            return False, "ZERO_LENGTH_SEQUENCE"
        if record.id in sequence_names:
            return False, "DUPLICATE_SEQUENCE_NAME"
        sequence_names.add(record.id)
    return True, "VALID"


def warn_parasail(filepath):
    try:
        records = list(SeqIO.parse(filepath, "fasta"))
    except Exception:
        return False

    if len(records) > 500:
        return True

    for record in records:
        if len(record.seq) > 20000:
            return True

    return False
