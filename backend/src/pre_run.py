import re
from Bio.SeqUtils import gc_fraction


def run_preprocessing(raw_seq_dict):
    seq_dict = {}
    for key, record in raw_seq_dict.items():
        sequence = str(record.seq).upper()
        sequence = re.sub(r"[^A-Z]", "", sequence)
        key= re.sub(r"[, -]", "_", key)
        key = key[:50] ## hardcap for id length
        seq_dict[key] = str(sequence)
    return seq_dict

def grab_stats(seq_dict):
    seq_stats = {}
    for key, record_str in seq_dict.items():
        gcCount = round(gc_fraction(record_str), 2) * 100
        genLen = len(record_str)
        seq_stats.setdefault(key, [])
        seq_stats[key].append(gcCount)
        seq_stats[key].append(genLen)
    return seq_stats

def residue_check(seq):
    return bool(re.search(r"[EFILPQZ]", seq))