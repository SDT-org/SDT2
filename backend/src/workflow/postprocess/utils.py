from typing import Iterable, Dict
from Bio.SeqRecord import SeqRecord
from Bio.SeqUtils import gc_fraction


def get_seq_stats(seq_dict: Dict[str, SeqRecord]) -> dict[str, list[float]]:
    seq_stats = {}
    for id, seq in seq_dict.items():
        gcCount = round(gc_fraction(seq), 2) * 100
        genLen = len(str(seq))
        seq_stats.setdefault(id, [])
        seq_stats[id].append(gcCount)
        seq_stats[id].append(genLen)
    return seq_stats
