import subprocess
import os
import pandas as pd
import numpy as np
from heatmap import dataframe_to_triangle
from pre_run import run_preprocessing, grab_stats, residue_check, preprocess_fasta_to_fasta
from document_paths import build_document_paths
from process_data import save_stats_to_csv, seq_dict_to_json, output_summary, save_matrix_to_csv
import time
from datetime import datetime
from Bio import SeqIO
import shutil
from pandas import DataFrame

def lzani_to_full_matrix(results_tsv_path, ids_tsv_path, score_column='ani'):
    all_ids = pd.read_csv(ids_tsv_path, sep='\t')['id'].tolist()
    df = pd.read_csv(results_tsv_path, sep='\t')
    matrix = df.pivot_table(index='query', columns='reference',
                           values=score_column, aggfunc='first') * 100
    matrix = matrix.reindex(index=all_ids, columns=all_ids) ## can probs plug in reorder logic here
    matrix = matrix.fillna(25.0) ## filling NAs with biological floor
    matrix = matrix.replace(0.0, 25.0)
    symmetric_matrix = (matrix + matrix.T) / 2 ## average the two triangles for symmetic matrix??
    np.fill_diagonal(symmetric_matrix.values, 100.0)
    return symmetric_matrix

def get_output_paths(output_prefix_base, out_dir):
    os.makedirs(out_dir, exist_ok=True)
    prefix = os.path.basename(output_prefix_base or "lzani_out")
    return (
        os.path.join(out_dir, f"{prefix}_results.tsv"),
        os.path.join(out_dir, f"{prefix}_ids.tsv")
    )

def run_lzani(
    settings: dict,
    set_stage=lambda stage: print(f"Stage: {stage}")
):
    def lzani_all2all(lz_ani_executable, input_fasta,
                      results_tsv_path, ids_tsv_path, threads):
        abs_lzani_executable_path = os.path.abspath(os.path.normpath(lz_ani_executable))
        cmd = [
            abs_lzani_executable_path, "all2all",
            "--in-fasta", input_fasta,
            "--out", results_tsv_path, "--out-ids", ids_tsv_path,
            "--out-format", "complete", "--threads", str(threads),
            "--verbose", "2"
        ]
        process = subprocess.run(cmd, check=True, capture_output=True, text=True)

    start_time, start_counter = datetime.now(), time.perf_counter()
    input_file = settings["input_file"]
    file_name = os.path.basename(input_file)
    out_dir = settings["out_dir"]
    doc_paths = build_document_paths(out_dir)
    set_stage("Preprocessing")
    preprocessed_fasta_path = os.path.join(out_dir, "preprocessed.fasta")
    preprocess_fasta_to_fasta(input_file, preprocessed_fasta_path)

    sequences = SeqIO.parse(open(preprocessed_fasta_path, encoding="utf-8"), "fasta")
    raw_seq_dict_for_preprocessing = {record.id: record for record in sequences}
    seq_dict_processed_strings = run_preprocessing(raw_seq_dict_for_preprocessing)
    seq_stats = grab_stats(seq_dict_processed_strings)

    set_stage("Analyzing with LZ-ANI")
    results_tsv_path, ids_tsv_path = get_output_paths(
        settings.get("lz_ani_output_prefix"),
        out_dir
    )
    lzani_all2all(
        lz_ani_executable=settings["lz_ani_executable_path"],
        input_fasta=preprocessed_fasta_path,
        results_tsv_path=results_tsv_path,
        ids_tsv_path=ids_tsv_path,
        threads=settings["num_processes"]
    )
    matrix = lzani_to_full_matrix(
        results_tsv_path,
        ids_tsv_path,
        score_column=settings["score_type"]
    )

    cluster_method = settings.get("cluster_method")
    if cluster_method and cluster_method != "None":
        set_stage("Clustering")
        from scipy.cluster import hierarchy
        dist_matrix = 100.0 - matrix.values
        np.fill_diagonal(dist_matrix, 0)
        condensed_dist = [dist_matrix[i, j] for i in range(len(matrix)) for j in range(i + 1, len(matrix))]
      ##--need to plug in the reorder logic from parasail here as a component
        linked = hierarchy.linkage(np.array(condensed_dist), method=cluster_method)
        dendro_data = hierarchy.dendrogram(linked, orientation='right', no_plot=True)
        reordered_indices = [int(i) for i in dendro_data['ivl']]
        seq_ids = matrix.index.tolist()
        reordered_ids = [seq_ids[i] for i in reordered_indices]
        matrix = matrix.loc[reordered_ids, reordered_ids]
            
    set_stage("Postprocessing")
    results_df = pd.read_csv(results_tsv_path, sep='\t')
    results_df.to_csv(doc_paths.columns, sep=',', index=False)
    save_stats_to_csv(seq_stats, doc_paths.stats)
    df = DataFrame(matrix, index=matrix.index, columns=matrix.index)
    save_matrix_to_csv(df, doc_paths.matrix, doc_paths.triangle)
    seq_dict_to_json(seq_dict_processed_strings, doc_paths.seq_dict)

    set_stage("Finalizing")
    end_time, end_counter = datetime.now(), time.perf_counter()
    summary_text = output_summary(file_name, start_time, end_time, start_counter, end_counter)
    with open(doc_paths.summary, "w") as file: file.write(summary_text)
