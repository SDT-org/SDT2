from typing import Iterable
from Bio.SeqRecord import SeqRecord
from Bio.SeqUtils import gc_fraction


def get_seq_stats(records: Iterable[SeqRecord]):
    seq_stats = {}
    for record in records:
        gcCount = round(gc_fraction(record.seq), 2) * 100
        genLen = len(str(record.seq))
        seq_stats.setdefault(record.id, [])
        seq_stats[record.id].append(gcCount)
        seq_stats[record.id].append(genLen)
    return seq_stats
