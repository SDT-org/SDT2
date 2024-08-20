import time
import platform
import psutil
import os
import sys
import re
import random
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from itertools import combinations_with_replacement as cwr
from functools import partial
import numpy as np
import pandas as pd
import parasail
from Bio import SeqIO, Phylo
from Bio.Align import PairwiseAligner
from Bio.Phylo.TreeConstruction import (
    DistanceTreeConstructor,
    DistanceMatrix,
)

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "../../")))
from config import app_version

start_time = time.time()
start_run = time.strftime("%Y-%m-%d %H:%M:%S", time.localtime())


def make_aligner(settings):
    is_aa = settings["is_aa"]
    aligner = PairwiseAligner()
    aligner.match_score = getattr(settings, "match", 1.5)
    aligner.mismatch_score = getattr(settings, "mismatch", -1)
    aligner.target_internal_open_gap_score = getattr(
        settings, "iog", -2 if is_aa else -3
    )
    aligner.target_internal_extend_gap_score = getattr(settings, "ieg", -2)
    aligner.target_left_open_gap_score = getattr(settings, "log", -3 if is_aa else -4)
    aligner.target_left_extend_gap_score = getattr(settings, "leg", -2 if is_aa else -3)
    aligner.target_right_open_gap_score = getattr(settings, "rog", -3 if is_aa else -4)
    aligner.target_right_extend_gap_score = getattr(
        settings, "reg", -2 if is_aa else -3
    )
    aligner.query_internal_open_gap_score = getattr(
        settings, "iog", -2 if is_aa else -3
    )
    aligner.query_internal_extend_gap_score = getattr(settings, "ieg", -2)
    aligner.query_left_open_gap_score = getattr(settings, "log", -3 if is_aa else -4)
    aligner.query_left_extend_gap_score = getattr(settings, "leg", -2 if is_aa else -3)
    aligner.query_right_open_gap_score = getattr(settings, "rog", -3 if is_aa else -4)
    aligner.mode = settings["alignment_type"]
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


def process_pair(id_sequence_pair, settings):
    ids = id_sequence_pair[0]
    id1 = ids[0]
    id2 = ids[1]
    seq1 = id_sequence_pair[1][0]
    seq2 = id_sequence_pair[1][1]
    aligner = make_aligner(settings)

    aln = parasail.nw_trace(seq1, seq2, 10, 1, parasail.blosum62)
    score = get_similarity(aln.traceback.query, aln.traceback.ref)

    if settings.get("aln_out"):
        fname = os.path.join(
            settings["aln_out"], str(id1 + "__" + id2) + "_aligned.fasta"
        )
        aligned_sequences = list(aln)
        seq_records = [
            SeqRecord(Seq(aligned_sequences[i]), id=ids[i], description="")
            for i in range(len(aligned_sequences))
        ]
        SeqIO.write(seq_records, fname, "fasta")

    return id_sequence_pair[0], score


## Creates arrays for storing scores as matrix and set pool for multiprocessing.
def get_alignment_scores(
    seq_dict, settings, pool, cancelled, increment_pair_progress, set_pair_count
):
    seq_ids = list(seq_dict.keys())
    n = len(seq_ids)
    dist_scores = np.zeros((n, n))
    # for each sequence id in seq_id add to new dict as key and the index position is enumerated and stored as value for later reference
    order = {seq_id: i for i, seq_id in enumerate(seq_ids)}
    # create list  combinations including self v. self
    combos = list(cwr(seq_ids, 2))
    id_sequence_pairs = []

    for ids in combos:
        id_sequence_pairs.append([ids, [seq_dict[ids[0]], seq_dict[ids[1]]]])

    # calculate the total combinations witout self
    total_pairs = sum(1 for _ in cwr(seq_ids, 2))
    set_pair_count(total_pairs)

    print(f"\rNumber of sequences: {len(seq_ids)}\r", flush=True)
    print(f"\rNumber of pairs: {total_pairs}\r", flush=True)

    # Randomly sample if sequence contains amino acid residues
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


# Reorder the scores matrix based on the tree and save it to a new CSV
def tree_clustering(args, dm, filename):
    constructor = DistanceTreeConstructor()
    clustering_method = getattr(constructor, args["cluster_method"])
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


def output_summary(args, file_name, start_time, end_time, run_summary):
    aligner = make_aligner(args)
    aligner_params_str = str(aligner).split("\n")
    # Adjust the indentation of each line
    formatted_aligner_params = "\n    ".join(aligner_params_str)
    build_type = f"{platform.system()} {platform.release()} {platform.machine()}"
    kernel_info = platform.processor()
    cpu_cores = os.cpu_count()
    total_ram = psutil.virtual_memory().total / (1024**3)

    # Run Summary

    output_content = f"""
    SDT {app_version} release for {platform.system()}
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


def process_data(
    settings,
    pool,
    cancelled,
    increment_pair_progress=lambda: None,
    set_stage=lambda stage: print(f"Stage: {stage}"),
    set_pair_count=lambda count: print(f"Pair count: {count}"),
):
    input_file = settings["input_file"]
    INPUT_PATH = input_file

    file_name = os.path.basename(input_file)
    file_base = os.path.splitext(file_name)[0]

    set_stage("Preparing")
    print("Stage: Preparing")

    sequences = SeqIO.parse(open(INPUT_PATH, encoding="utf-8"), "fasta")

    # Create a dictionary with the full description as the key, replacing spaces with underscores
    processed_seq_dict = {
        record.description.replace(" ", ""): record for record in sequences
    }

    set_stage("Preprocessing")
    seq_dict = run_preprocessing(processed_seq_dict)

    set_stage("Analyzing")
    print("Stage: Analyzing")

    dist_scores, order = get_alignment_scores(
        seq_dict, settings, pool, cancelled, increment_pair_progress, set_pair_count
    )

    if cancelled.value:
        return

    set_stage("Postprocessing")
    print("Stage: Postprocessing")

    order = list(order.keys())

    dm = create_distance_matrix(dist_scores, order)
    aln_scores = 100 - dist_scores

    out_dir = settings["out_dir"]
    cluster_method = settings["cluster_method"]

    if cluster_method == "nj" or cluster_method == "upgma":
        _, new_order = tree_clustering(
            settings, dm, os.path.join(out_dir, f"{file_base}")
        )
        reorder_index = [
            order.index(id_) for id_ in new_order
        ]  # create numerical index of order and  of new order IDs

        aln_reordered = aln_scores[reorder_index, :][:, reorder_index]
        aln_lowt = np.tril(np.around(aln_reordered, 2))
        aln_lowt[np.triu_indices(aln_lowt.shape[0], k=1)] = np.nan
        # Create a DataFrame from the lower triangular matrix
        df = pd.DataFrame(aln_lowt, index=new_order)
        save_matrix_to_csv(df, os.path.join(out_dir, f"{file_base}_mat.csv"))
        # save_cols_to_csv(df, new_order, os.path.join(args.out_dir, f"{file_base}"))
    else:
        aln_lowt = np.tril(np.around(aln_scores, 2))
        aln_lowt[np.triu_indices(aln_lowt.shape[0], k=1)] = np.nan
        df = pd.DataFrame(aln_lowt, index=order)
        save_matrix_to_csv(df, os.path.join(out_dir, f"{file_base}_mat.csv"))
        # save_cols_to_csv(df, new_order, os.path.join(args.out_dir, f"{file_base}"))

    set_stage("Finalizing")
    print("Stage: Finalizing")

    # Finalize run
    end_time = time.time()
    end_run = time.strftime("%Y-%m-%d %H:%M:%S", time.localtime())
    elapsed_time = end_time - start_time
    run_summary = time.strftime("%H:%M:%S", time.gmtime(elapsed_time))
    save_output_summary = output_summary(
        settings,
        file_name,
        start_run,
        end_run,
        run_summary,
    )
    # Write to a text file
    with open(os.path.join(out_dir, f"{file_base}_summary.txt"), "w") as file:
        file.write(save_output_summary)
    print(f"Elapsed time: {elapsed_time} seconds")
