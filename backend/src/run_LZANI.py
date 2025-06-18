import subprocess
import os
import pandas as pd
import numpy as np
from heatmap import dataframe_to_triangle

def lzani_to_full_matrix(results_tsv_path, ids_tsv_path, score_column='ani'):
    all_ids = pd.read_csv(ids_tsv_path, sep='\t')['id'].tolist()
    
    df = pd.read_csv(results_tsv_path, sep='\t')
    matrix = df.pivot_table(index='query', columns='reference', 
                           values=score_column, aggfunc='first') * 100
    
    matrix = matrix.reindex(index=all_ids, columns=all_ids)
    
    matrix = matrix.fillna(25.0)
    matrix = matrix.replace(0.0, 25.0)
    symmetric_matrix = (matrix + matrix.T) / 2
    
    np.fill_diagonal(symmetric_matrix.values, 100.0)
    
    return symmetric_matrix

def get_output_paths(output_prefix_base, out_dir):
    os.makedirs(out_dir, exist_ok=True)
    prefix = os.path.basename(output_prefix_base or "lzani_out")
    return (
        os.path.join(out_dir, f"{prefix}_results.tsv"),
        os.path.join(out_dir, f"{prefix}_ids.tsv")
    )
    
def execute_lzani(lz_ani_executable, input_fasta, 
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
    
def run_lzani(
    raw_input_fasta: str,
    lz_ani_executable_path: str,
    score_type: str,
    out_dir: str,
    threads: int = 0,
    lz_ani_output_prefix: str | None = None
):
    results_tsv_path, ids_tsv_path = get_output_paths(lz_ani_output_prefix, out_dir)
    
    execute_lzani(
        lz_ani_executable=lz_ani_executable_path,
        input_fasta=raw_input_fasta,
        results_tsv_path=results_tsv_path,
        ids_tsv_path=ids_tsv_path,
        threads=threads
    )
    
    matrix = lzani_to_full_matrix(results_tsv_path, ids_tsv_path, score_column=score_type)
    
    return matrix
