import mimetypes
import re
from typing import Iterator
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from workflow.models import WorkflowResult

INVALID_SYMBOLS_RE = re.compile(r"[, -]")
INVALID_CHARACTERS_RE = re.compile(r"[^A-Z]")
RESIDUE_RE = re.compile(r"[EFILPQZ]")


def run(result: WorkflowResult, fasta_path: str) -> WorkflowResult:
    filetype, _ = mimetypes.guess_type(fasta_path)

    if not filetype or filetype != "text/fasta":
        result.errors.append("INVALID_FILE_TYPE")
        return result

    try:
        unparsed_records = SeqIO.parse(fasta_path, "fasta")
    except Exception as exception:
        result.errors.append(f"INVALID_FASTA_FORMAT: {exception}")
        return result

    result = parse_records(unparsed_records, result)

    return result


def format_sequence_record(record: SeqRecord) -> SeqRecord:
    formatted_id = re.sub(r"[, -]", "_", record.id or "")[:50]
    upcased_sequence = str(record.seq).upper()
    formatted_sequence = re.sub(r"[^A-Z]", "", upcased_sequence)
    return SeqRecord(Seq(formatted_sequence), id=formatted_id, description="")


def residue_check(seq):
    return bool(RESIDUE_RE.search(seq))


def parse_records(input: Iterator[SeqRecord], result: WorkflowResult) -> WorkflowResult:
    if result.seq_dict:
        # If records are already parsed, we can skip parsing again
        return result

    max_length = 0
    count = 0
    is_aa = result.is_aa
    sequence_names = set()
    seq_dict = result.seq_dict or {}

    for record in input:
        count += 1

        if not record.id:
            result.errors.append("MISSING_SEQUENCE_ID")
            return result

        if not record.seq:
            result.errors.append("ZERO_LENGTH_SEQUENCE")
            return result

        if record.id in sequence_names:
            result.errors.append("DUPLICATE_SEQUENCE_NAME")
            return result

        sequence_names.add(record.id)

        # replace some characters so we can export to CSV later
        formatted_id = INVALID_SYMBOLS_RE.sub("_", record.id)[:50]
        upcased_sequence = str(record.seq).upper()
        formatted_sequence = INVALID_CHARACTERS_RE.sub("", upcased_sequence)

        if not formatted_sequence:
            result.errors.append("ZERO_LENGTH_SEQUENCE")
            return result

        # check every three sequences for amino acids
        # TODO: do we want to check this every time?
        if is_aa == None and count % 3 == 0:
            is_aa = residue_check(formatted_sequence)

        max_length = max(max_length, len(record.seq))
        seq_dict[formatted_id] = formatted_sequence

    if count < 2:
        result.errors.append("NOT_ENOUGH_SEQUENCES")

    return result._replace(
        seq_dict=seq_dict,
        max_sequence_length=max_length,
        is_aa=is_aa,
    )
