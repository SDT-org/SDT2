import time
import platform
import psutil
import os
import sys
import random
import json
import numpy as np
from datetime import datetime
from pandas import DataFrame
from numpy import zeros
from Bio import SeqIO
from scipy.cluster import hierarchy
from heatmap import dataframe_to_triangle
from document_paths import build_document_paths
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "../../")))
from config import app_version
from pre_run import run_preprocessing, grab_stats
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "analysis")))
from run_parasail import get_alignment_scores




def save_matrix_to_csv(df, matrix_path, triangle_path):
    index = df.index
    try:
        triangle = dataframe_to_triangle(df)
        triangle.index = index
        triangle.to_csv(triangle_path, mode="wt", header=False, index=True, sep=',')
    except Exception:
        pass
    df.to_csv(matrix_path, mode="w", header=False, index=True, sep=',')

def save_cols_to_csv(df, path):
    order = df.index
    df.columns = df.index
    columnar_output = []
    for i_idx, row_id in enumerate(order):
        for j_idx, col_id in enumerate(order):
            if i_idx > j_idx:
                columnar_output.append([row_id, col_id, df.loc[row_id, col_id]])
    columnar_df = DataFrame(columnar_output, columns=["First Sequence", "Second Sequence", "Identity Score"])
    columnar_df.to_csv(path, mode="w", header=True, index=False, sep=',')

def save_stats_to_csv(seq_stats, filename):
    stats_list = []
    for key, value in seq_stats.items():
        stats_list.append([key, value[0], value[1]])
    stats_df = DataFrame(stats_list, columns=["Sequence", "GC %", "Sequence Length"])
    stats_df.to_csv(filename, mode="w", header=True, index=False, sep=',')

def seq_dict_to_json(seq_dict, filename):
    with open(filename,"w") as file: json.dump(seq_dict, file, indent=4)

def friendly_total_time(total_time):
    m, s = divmod(total_time, 60)
    return f'{int(m)}m {s:.2f}s' if m > 0 else f'{s:.2f}s'

def output_summary(file_name, start_time, end_time, start_counter, end_counter):
    build_type = f"{platform.system()} {platform.release()} {platform.machine()}"
    total_cores, total_ram = os.cpu_count(), psutil.virtual_memory().total / (1024**3)
    total_counter = end_counter - start_counter
    return f"""
    SDT {app_version} release for {platform.system()}
    Developed by Michael Lund and Josiah Ivey; Brejnev Muhire,
    Darren Martin, Simona Kraberger, Qiyun Zhu, Pierre Lefeuvre, Jean-Michele Lett, Philippe Roumagnac, Arvind Varsani
    System info: Host: {platform.node()}, OS: {build_type}, CPU: {platform.processor()} - {total_cores} cores, Memory: {total_ram:.2f} GB RAM
    Run info for {file_name}: Start: {start_time.strftime("%b %d %Y, %I:%M %p %Z")}, End: {end_time.strftime("%b %d %Y, %I:%M %p %Z")}, Total: {friendly_total_time(total_counter)}
    Parasail: Using Needleman-Wunsch (stats) algorithm. Nucleotide: Open={13}, Extend={1} (BLOSUM62). Amino acid: Open={10}, Extend={1} (BLOSUM62).
    """

def fasta_alignments(seq_records, fname):
    SeqIO.write(seq_records, fname, "fasta")

def process_data(
    settings, pool, cancelled,
    increment_pair_progress=lambda: None,
    set_stage=lambda stage: print(f"Stage: {stage}"),
    set_pair_count=lambda count: print(f"Pair count: {count}"),
):
    start_time, start_counter = datetime.now(), time.perf_counter()
    input_file, file_name = settings["input_file"], os.path.basename(settings["input_file"])

    set_stage("Preparing")
    sequences = SeqIO.parse(open(input_file, encoding="utf-8"), "fasta")
    raw_seq_dict_for_preprocessing = {record.id: record for record in sequences}

    set_stage("Preprocessing")

    seq_dict_processed_strings = run_preprocessing(raw_seq_dict_for_preprocessing)
    seq_stats = grab_stats(seq_dict_processed_strings)

    set_stage("Analyzing")
    dist_scores, order_map = get_alignment_scores(
        seq_dict_processed_strings, settings, pool, cancelled, increment_pair_progress, set_pair_count
    )

    if cancelled.value: return

    set_stage("Clustering")
    seq_ids_in_order = list(order_map.keys())

    aln_scores = 100.0 - dist_scores

    out_dir, cluster_method = settings["out_dir"], settings.get("cluster_method")
    doc_paths = build_document_paths(out_dir)

    final_ordered_ids, final_matrix_for_df = seq_ids_in_order, aln_scores
##3 lets bring this out as a compoennt for LZAZANI and parasail alike, just do some handling
    if cluster_method is not None and len(seq_ids_in_order) > 1:
        current_dist_matrix_for_clustering = np.copy(dist_scores)
        np.fill_diagonal(current_dist_matrix_for_clustering, 0)
        condensed_dist_matrix = [current_dist_matrix_for_clustering[i, j] for i in range(len(seq_ids_in_order)) for j in range(i + 1, len(seq_ids_in_order))]

        if condensed_dist_matrix:
            try:
                linked = hierarchy.linkage(np.array(condensed_dist_matrix), method=cluster_method)
                dendro_data = hierarchy.dendrogram(linked, orientation='right', no_plot=True)
                reordered_original_indices = [int(i_str) for i_str in dendro_data['ivl']]
                final_ordered_ids = [seq_ids_in_order[i] for i in reordered_original_indices]
                final_matrix_for_df = aln_scores[np.ix_(reordered_original_indices, reordered_original_indices)]
            except Exception as e:
                print(f"Warning: Clustering failed: {e}. Using original order.", file=sys.stderr)
                final_ordered_ids, final_matrix_for_df = seq_ids_in_order, aln_scores
        else:
            final_ordered_ids, final_matrix_for_df = seq_ids_in_order, aln_scores
##
    df = DataFrame(final_matrix_for_df, index=final_ordered_ids, columns=final_ordered_ids)

    set_stage("Postprocessing")
    save_cols_to_csv(df, doc_paths.columns)
    save_stats_to_csv(seq_stats, doc_paths.stats)
    save_matrix_to_csv(df, doc_paths.matrix, doc_paths.triangle)
    seq_dict_to_json(seq_dict_processed_strings, doc_paths.seq_dict)

    set_stage("Finalizing")
    end_time, end_counter = datetime.now(), time.perf_counter()
    summary_text = output_summary(file_name, start_time, end_time, start_counter, end_counter)
    with open(doc_paths.summary, "w") as file: file.write(summary_text)
    print(f"Elapsed time: {friendly_total_time(end_counter - start_counter)}")
