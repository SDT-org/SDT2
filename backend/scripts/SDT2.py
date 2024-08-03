import time
import datetime
import platform
import psutil
import os

# import logging
import sys
import re
import random
import argparse
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from itertools import combinations_with_replacement as cwr
from multiprocessing import Pool, Manager, cpu_count
from functools import partial
import numpy as np
import pandas as pd
from Bio import SeqIO, Phylo
from Bio.Align import PairwiseAligner
from Bio.Phylo.TreeConstruction import (
    DistanceTreeConstructor,
    DistanceMatrix,
)

start_time = time.time()
print("time start")
start_run = time.strftime("%Y-%m-%d %H:%M:%S", time.localtime())
useable_processes = cpu_count() - 1
aligner_params = None
aligner_alg = None


# Arguments for run
def parse_args():
    parser = argparse.ArgumentParser(description="Argument parser for SDT2")
    parser.add_argument(
        "--cluster_method",
        type=str,
        choices=["nj", "upgma", "None"],
        help="order output by clusting algorithm",
    )
    parser.add_argument(
        "--out_dir",
        type=str,
        default=os.path.dirname(os.path.abspath(__file__)),
        help="the output directory",
    )
    parser.add_argument(
        "--alignment_type",
        type=str,
        default="global",
        help="alignment type global/local",
    )
    parser.add_argument("input_file", type=str, help="The input fasta file")
    parser.add_argument(
        "--num_processes",
        type=int,
        default=useable_processes,
        help="number of processes to run",
    )
    parser.add_argument(
        "--aln_out", type=str, default=None, help="Output directory for"
    )
    parser.add_argument("--match", type=float, help="aligner match score")
    parser.add_argument("--mismatch", type=int, help="aligner mismatch score")
    parser.add_argument("--iog", type=int, help="internal open gap score")
    parser.add_argument("--ieg", type=int, help="internal left extend gap score")
    parser.add_argument("--log", type=int, help="left open gap score")
    parser.add_argument("--leg", type=int, help="left extend gap score")
    parser.add_argument("--rog", type=int, help="right open gap score")
    parser.add_argument("--reg", type=int, help="right extend gap score")
    return parser.parse_args()

def make_aligner(is_aa = False):
    args = parse_args()
    aligner = PairwiseAligner()
    aligner.match_score = getattr(args, "match", 1.5)
    aligner.mismatch_score = getattr(args, "mismatch", -1)
    aligner.target_internal_open_gap_score = getattr(args, "iog", -2 if is_aa else -3)
    aligner.target_internal_extend_gap_score = getattr(args, "ieg", -2)
    aligner.target_left_open_gap_score = getattr(args, "log", -3 if is_aa else -4)
    aligner.target_left_extend_gap_score = getattr(args, "leg", -2 if is_aa else -3)
    aligner.target_right_open_gap_score = getattr(args, "rog", -3 if is_aa else -4)
    aligner.target_right_extend_gap_score = getattr(args, "reg", -2 if is_aa else -3)
    aligner.query_internal_open_gap_score = getattr(args, "iog", -2 if is_aa else -3)
    aligner.query_internal_extend_gap_score = getattr(args, "ieg", -2)
    aligner.query_left_open_gap_score = getattr(args, "log", -3 if is_aa else -4)
    aligner.query_left_extend_gap_score = getattr(args, "leg", -2 if is_aa else -3)
    aligner.query_right_open_gap_score = getattr(args, "rog", -3 if is_aa else -4)
    aligner.query_right_extend_gap_score = getattr(args, "reg", -2 if is_aa else -3)
    aligner.mode = args.alignment_type
    return aligner

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


# Performs sequence alignment using BioPython PairwiseAligner
def biopython_align(id_sequence_pair, is_aa):
    args = parse_args()
    ids = id_sequence_pair[0]
    id1 = ids[0]
    id2 = ids[1]
    seq1 = id_sequence_pair[1][0]
    seq2 = id_sequence_pair[1][1]
    aligner = make_aligner(is_aa)

    aln = aligner.align(seq1, seq2)[0]
    score = get_similarity(aln[0], aln[1])

    if args.aln_out:
        fname = os.path.join(args.aln_out, str(id1 + "__" + id2) + "_aligned.fasta")
        aligned_sequences = list(aln)
        seq_records = [
            SeqRecord(Seq(aligned_sequences[i]), id=ids[i], description="")
            for i in range(len(aligned_sequences))
        ]  ## rewrite for multi
        SeqIO.write(seq_records, fname, "fasta")

    return score


def process_pair(id_sequence_pair, counter, total_pairs, is_aa):
    score = biopython_align(id_sequence_pair, is_aa)
    counter.value += 1
    print(
        "\rPerforming alignment: progress "
        + str(int((float(counter.value) / total_pairs) * 100))
        + "% - pair",
        +counter.value,
        flush=True,
    )
    return id_sequence_pair[0], score


## Creates arrays for storing scores as matrix and set pool for multiprocessing.
def get_alignment_scores(seq_dict, args):
    num_processes = args.num_processes
    seq_ids = list(seq_dict.keys())
    n = len(seq_ids)
    dist_scores = np.zeros((n, n))
    # for each sequence id in seq_id add to new dict as key and the index position is enumerated and stored as value for later reference
    order = {seq_id: i for i, seq_id in enumerate(seq_ids)}
    manager = Manager()
    counter = manager.Value("i", 0)
    # create list  combinations including self v. self
    combos = list(cwr(seq_ids, 2))
    id_sequence_pairs = []

    for ids in combos:
        id_sequence_pairs.append([ids, [seq_dict[ids[0]], seq_dict[ids[1]]]])

    # calculate the total combinations witout self
    total_pairs = sum(1 for _ in cwr(seq_ids, 2))
    print(f"\rNumber of sequences: {len(seq_ids)}\r", flush=True)
    print(f"\rNumber of pairs: {total_pairs}\r", flush=True)

    # Randomly sample if sequence contains amino acid residues
    sample_size = 3
    sampled_seqs = random.sample(
        list(seq_dict.values()), min(len(seq_dict), sample_size)
    )
    is_aa = any(residue_check(seq) for seq in sampled_seqs)

    bound_process_pair = partial(
        process_pair,
        counter=counter,
        total_pairs=total_pairs,
        is_aa=is_aa,
    )

    with Pool(processes=num_processes) as pool:
        results = pool.map(bound_process_pair, id_sequence_pairs)
        for result in results:
            [seqid1, seqid2], score = result
            seqid1, seqid2 = order[seqid1], order[seqid2]
            dist_scores[seqid1, seqid2] = score
            dist_scores[seqid2, seqid1] = score
    return dist_scores, order


# Reorder the scores matrix based on the tree and save it to a new CSV
def tree_clustering(dm, filename):
    args = parse_args()
    constructor = DistanceTreeConstructor()
    clustering_method = getattr(constructor, args.cluster_method)
    tree_file = filename + "_tree.nwk"
    tree = clustering_method(dm)
    Phylo.write(tree, tree_file, "newick")
    tree = Phylo.read(tree_file, "newick")
    new_order = [leaf.name for leaf in tree.get_terminals()]
    return tree_file, new_order  # Return the tree file name instead of the tree object


# Format distance scores to triangle matrix for NJ tree creation
def create_distance_matrix(dist_scores, order):
    lower_triangle = [
        dist_scores[i, : i + 1].tolist() for i in range(len(order))
    ]  ##this has to be list, cant use np
    # Create the DistanceMatrix object
    dm = DistanceMatrix(order, lower_triangle)
    return dm


# Save similarity scores as 2d matrix csv
def save_matrix_to_csv(df, filename):
    df.to_csv(filename, mode="w", header=False, index=True)


# Save similarity scores as 3 column csv
# def save_cols_to_csv(df, order, filename):
#     columnar_output = []
#     for i, row in enumerate(order):
#         for j, col in enumerate(order):
#             if i > j:  # lower triangular part (excluding diagonal)
#                 columnar_output.append([row, col, df.loc[row, col]])

#     # Convert to a DataFrame
#     columnar_df = pd.DataFrame(
#         columnar_output,
#         columns=["First Sequence", "Second Sequence", "Identity Score"],
#     )
#     columnar_df.to_csv(filename + "_cols.csv", mode="w", header=True, index=False)


def output_summary(file_name, start_time, end_time, run_summary):
    aligner = make_aligner()
    aligner_params_str = str(aligner).split("\n")
    # Adjust the indentation of each line
    formatted_aligner_params = "\n    ".join(aligner_params_str)
    build_type = f"{platform.system()} {platform.release()} {platform.machine()}"
    kernel_info = platform.processor()
    cpu_cores = os.cpu_count()
    total_ram = psutil.virtual_memory().total / (1024**3)

    # Run Summary

    output_content = f"""
    SDT v2.0.0-beta release for {platform.system()} built 14/06/2024
    Developed by Michael Lund, Josiah Ivey, Brejnev Muhire,
    Darren Martin, Simona Kraberger, Qiyun Zhu, Pierre Lefeuvre, Philippe Roumagnac, Arvind Varsani

    System info:
    Host:    {platform.node()} ({kernel_info}, {total_ram:.2f} GB RAM)
    Kernel:  {kernel_info} - auto-detect threads ({cpu_cores} CPU cores detected)
    OS - build: {build_type}
    Run info for {file_name}:
    Start time:    {start_time}

    BioPython {formatted_aligner_params}
    Using the {str(aligner.algorithm)}

    End time: {end_time}
    Total runtime: {run_summary}
    """

    return output_content


def fasta_alignments(seq_records, fname):
    SeqIO.write(seq_records, fname, "fasta")


def main():
    args = parse_args()
    INPUT_PATH = args.input_file
    # logging.basicConfig(level=logging.DEBUG)

    file_path = INPUT_PATH
    file_name = os.path.basename(args.input_file)
    file_base = os.path.splitext(file_name)[0]

    print("Stage: Preparing")

    sequences = SeqIO.parse(open(INPUT_PATH, encoding="utf-8"), "fasta")

    # Create a dictionary with the full description as the key, replacing spaces with underscores
    processed_seq_dict = {
        record.description.replace(" ", ""): record for record in sequences
    }

    seq_dict = run_preprocessing(processed_seq_dict)

    print("Stage: Analyzing")

    dist_scores, order = get_alignment_scores(seq_dict, args)

    order = list(order.keys())

    dm = create_distance_matrix(dist_scores, order)
    aln_scores = 100 - dist_scores

    if args.cluster_method == "nj" or args.cluster_method == "upgma":
        _, new_order = tree_clustering(dm, os.path.join(args.out_dir, f"{file_base}"))
        reorder_index = [
            order.index(id_) for id_ in new_order
        ]  # create numerical index of order and  of new order IDs

        aln_reordered = aln_scores[reorder_index, :][:, reorder_index]
        aln_lowt = np.tril(np.around(aln_reordered, 2))
        aln_lowt[np.triu_indices(aln_lowt.shape[0], k=1)] = np.nan
        # Create a DataFrame from the lower triangular matrix
        df = pd.DataFrame(aln_lowt, index=new_order)
        save_matrix_to_csv(df, os.path.join(args.out_dir, f"{file_base}_mat.csv"))
        # save_cols_to_csv(df, new_order, os.path.join(args.out_dir, f"{file_base}"))
    else:
        aln_lowt = np.tril(np.around(aln_scores, 2))
        aln_lowt[np.triu_indices(aln_lowt.shape[0], k=1)] = np.nan
        df = pd.DataFrame(aln_lowt, index=order)
        save_matrix_to_csv(df, os.path.join(args.out_dir, f"{file_base}_mat.csv"))
        # save_cols_to_csv(df, new_order, os.path.join(args.out_dir, f"{file_base}"))

    print("Stage: Finalizing")

    # Finalize run
    end_time = time.time()
    end_run = time.strftime("%Y-%m-%d %H:%M:%S", time.localtime())
    elapsed_time = end_time - start_time
    run_summary = time.strftime("%H:%M:%S", time.gmtime(elapsed_time))
    save_output_summary = output_summary(
        file_name,
        start_run,
        end_run,
        run_summary,
    )
    # Write to a text file
    with open(os.path.join(args.out_dir, f"{file_base}_summary.txt"), "w") as file:
        file.write(save_output_summary)
    print(f"Elapsed time: {elapsed_time} seconds")


if __name__ == "__main__":
    main()
