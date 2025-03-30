import time
import platform
import psutil
import os
import sys
import re
import random
from datetime import datetime
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from itertools import combinations_with_replacement as cwr
from functools import partial
from pandas import DataFrame
from numpy import tril, triu_indices, zeros, around, nan
import parasail
from Bio import SeqIO
from Bio.SeqUtils import gc_fraction
from cluster import get_linkage_method_order
from heatmap import dataframe_to_lower_triangle

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "../../")))
from config import app_version

def run_preprocessing(raw_seq_dict):
    # convert seqrecord element in raw dict to key,vale of ids and seqs
    seq_dict = {}
    for key, record in raw_seq_dict.items():
        # convert all seqs to upper
        sequence = str(record.seq).upper()
        # remove any non alphabetical characters
        sequence = re.sub(r"[^A-Z]", "", sequence)
        # TO DO -add exceptions for ambigious bases
        # add to new dict without seqrecord clutter

        seq_dict[key] = str(sequence)
    return seq_dict


def grab_stats(seq_dict):
    # https://stackoverflow.com/questions/20585920/how-to-add-multiple-values-to-a-dictionary-key
    seq_stats = {}
    for key, record in seq_dict.items():
        gcCount = round(gc_fraction(str(record)), 2) * 100
        genLen = len(record)
        seq_stats.setdefault(key, [])
        seq_stats[key].append(gcCount)
        seq_stats[key].append(genLen)
    return seq_stats


##Calculate the similarity scores from the alignments by iterating through each position in the alignemtn files as a zip
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


def residue_check(seq):
    return bool(re.search(r"[EFILPQZ]", seq))


def process_pair(id_sequence_pair, settings):
    id1 = id_sequence_pair[0][0]
    id2 = id_sequence_pair[0][1]
    seq1 = id_sequence_pair[1][0]
    seq2 = id_sequence_pair[1][1]
    if id1 == id2:
        return id_sequence_pair[0], 0
    open = 10 if settings["is_aa"] else 13

    try:
        result = parasail.nw_trace(seq1, seq2, open, 1, parasail.blosum62)
        query = result.traceback.query
    except:
        raise Exception("PARASAIL_TRACEBACK")

    score = get_similarity(query, result.traceback.ref)

    if settings.get("aln_out"):
        fname = os.path.join(
            settings["aln_out"], str(id1 + "_" + id2) + "_aligned.fasta"
        )
        seq_records = [
            SeqRecord(Seq(result.traceback.query), id=id1, description=""),
            SeqRecord(Seq(result.traceback.ref), id=id2, description=""),
        ]
        SeqIO.write(seq_records, fname, "fasta")

    return id_sequence_pair[0], score


## Creates arrays for storing scores as matrix and set pool for multiprocessing.
def get_alignment_scores(
    seq_dict, settings, pool, cancelled, increment_pair_progress, set_pair_count
):
    seq_ids = list(seq_dict.keys())
    n = len(seq_ids)
    dist_scores = zeros((n, n))
    # for each sequence id in seq_id add to new dict as key and the index position is enumerated and stored as value for later reference
    order = {seq_id: i for i, seq_id in enumerate(seq_ids)}
    # create list  combinations including self v. self
    combos = list(cwr(seq_ids, 2))
    id_sequence_pairs = []

    for ids in combos:
        id_sequence_pairs.append([ids, [seq_dict[ids[0]], seq_dict[ids[1]]]])

    # calculate the total combinations including self
    total_pairs = sum(1 for _ in cwr(seq_ids, 2))
    set_pair_count(total_pairs)

    print(f"\rNumber of sequences: {len(seq_ids)}\r", flush=True)
    print(f"\rNumber of pairs: {total_pairs}\r", flush=True)

    sample_size = 3
    sampled_seqs = random.sample(
        list(seq_dict.values()), min(len(seq_dict), sample_size)
    )
    is_aa = any(residue_check(seq) for seq in sampled_seqs)
    settings["is_aa"] = is_aa

    bound_process_pair = partial(
        process_pair,
        settings=settings,
    )

    with pool:
        results = pool.imap(bound_process_pair, id_sequence_pairs)
        for _, result in enumerate(results, 1):
            if cancelled.value:
                return dist_scores, order
            [seqid1, seqid2], score = result
            seqid1, seqid2 = order[seqid1], order[seqid2]
            dist_scores[seqid1, seqid2] = score
            dist_scores[seqid2, seqid1] = score
            increment_pair_progress()

    return dist_scores, order



# Save similarity scores as 2d matrix csv
def save_matrix_to_csv(df, filename):
    tri_matrix = dataframe_to_lower_triangle(df)
    tri_matrix.to_csv(filename + "_mat.csv", mode="w", header=False, index=True)
    df.to_csv("matrix.csv", mode="w", header=False, index=True)

# Save similarity scores as 3 column csv
def save_cols_to_csv(df, filename):
    order = df.index
    df.columns = df.index
    columnar_output = []
    for i, row in enumerate(order):
        for j, col in enumerate(order):
            if i > j:  # lower triangular part (excluding diagonal)
                columnar_output.append([row, col, df.loc[row, col]])
    # Convert to a DataFrame
    columnar_df = DataFrame(
        columnar_output,
        columns=["First Sequence", "Second Sequence", "Identity Score"],
    )
    columnar_df.to_csv(filename + "_cols.csv", mode="w", header=True, index=False)


def save_stats_to_csv(seq_stats, filename):
    stats_list = []
    for key, value in seq_stats.items():
        stats_list.append([key, value[0], value[1]])
    stats_df = DataFrame(stats_list, columns=["Sequence", "GC %", "Sequence Length"])
    stats_df.to_csv(filename + "_stats.csv", mode="w", header=True, index=False)

def friendly_total_time(total_time):
    m, s = divmod(total_time, 60)
    return f'{int(m)} minute, {s:.2f} second' if m == 1 else f'{int(m)} minutes, {s:.2f} seconds' if m > 0 else f'{s:.2f} seconds'

def output_summary(file_name, start_time, end_time, start_counter, end_counter):
    build_type = f"{platform.system()} {platform.release()} {platform.machine()}"
    total_cores = os.cpu_count()
    total_ram = psutil.virtual_memory().total / (1024**3)
    total_counter = end_counter - start_counter

    return f"""
    SDT {app_version} release for {platform.system()}
    Developed by Michael Lund and Josiah Ivey; Brejnev Muhire,
    Darren Martin, Simona Kraberger, Qiyun Zhu, Pierre Lefeuvre, Philippe Roumagnac, Arvind Varsani

    System info:
    Host:   {platform.node()}
    OS:     {build_type}
    CPU:    {platform.processor()} - {total_cores} cores
    Memory: {total_ram:.2f} GB RAM

    Run info for {file_name}:
    Start time: {start_time.strftime("%b %d %Y, %I:%M %p %Z")}
    End time:   {end_time.strftime("%b %d %Y, %I:%M %p %Z")}
    Total time: {friendly_total_time(total_counter)}

    Parasail
    Using the Needleman-Wunsch algorithm with affine gap scoring:
    Nucleotide: Open gap penalty: 13, Extend gap: 1
    Amino acid: Open gap penalty: 10, Extend gap: 1
    """

def fasta_alignments(seq_records, fname):
    SeqIO.write(seq_records, fname, "fasta")

def process_data(
    settings,
    pool,
    cancelled,
    increment_pair_progress=lambda: None,
    set_stage=lambda stage: print(f"Stage: {stage}"),
    set_pair_count=lambda count: print(f"Pair count: {count}"),
):
    start_time = datetime.now()
    start_counter = time.perf_counter()
    input_file = settings["input_file"]
    file_name = os.path.basename(input_file)
    file_base = os.path.splitext(file_name)[0]

    set_stage("Preparing")
    print("Stage: Preparing")

    sequences = SeqIO.parse(open(input_file, encoding="utf-8"), "fasta")

    # Create a dictionary with the full description as the key, replacing spaces with underscores
    processed_seq_dict = {
        record.description.replace(" ", ""): record for record in sequences
    }

    set_stage("Preprocessing")
    seq_dict = run_preprocessing(processed_seq_dict)
    seq_stats = grab_stats(seq_dict)

    set_stage("Analyzing")
    print("Stage: Analyzing")

    dist_scores, order = get_alignment_scores(
        seq_dict, settings, pool, cancelled, increment_pair_progress, set_pair_count
    )

    if cancelled.value:
        return

    set_stage("Clustering")
    print("Stage: Clustering")

    order = list(order.keys())

    aln_scores = 100 - dist_scores

    out_dir = settings["out_dir"]
    cluster_method = settings["cluster_method"]

    if cluster_method is not None:
        new_order = get_linkage_method_order(dist_scores, cluster_method, order)

        reorder_index = [
            order.index(id) for id in new_order
        ]  # create numerical index of order and  of new order IDs

        # numpy array indexing syntax rows/cols
        aln_scores = aln_scores[reorder_index, :][:, reorder_index]
        order = new_order

        df = DataFrame(aln_scores, index=new_order)

    df = DataFrame(aln_scores, index=order)

    set_stage("Postprocessing")
    print("Stage: Postprocessing")

    save_cols_to_csv(df, os.path.join(out_dir, file_base))
    save_stats_to_csv(seq_stats, os.path.join(out_dir, file_base))
    save_matrix_to_csv(df, os.path.join(out_dir, file_base))


    set_stage("Finalizing")
    print("Stage: Finalizing")

    end_time = datetime.now()
    end_counter = time.perf_counter()
    save_output_summary = output_summary(
        file_name,
        start_time,
        end_time,
        start_counter,
        end_counter
    )

    with open(os.path.join(out_dir, f"{file_base}_summary.txt"), "w") as file:
        file.write(save_output_summary)
    print(f"Elapsed time: {friendly_total_time(end_counter - start_counter)}")
