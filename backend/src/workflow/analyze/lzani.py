from multiprocessing.sharedctypes import Synchronized
import subprocess
import os
from typing import Callable
from workflow.models import RunSettings, WorkflowResult
from transformations import lzani_tsv_to_distance_matrix
import re


def run(
    result: WorkflowResult,
    settings: RunSettings,
    set_progress: Callable[[int], None],
    canceled: Synchronized,
) -> WorkflowResult:
    cmd = [
        os.path.abspath(os.path.normpath(settings.lzani.exec_path)),
        "all2all",
        "--in-fasta",
        settings.fasta_path,
        "--out",
        settings.doc_paths.lzani_results,
        "--out-ids",
        settings.doc_paths.lzani_results_ids,
        "--out-type", "tsv",
        "--out-format", "qidx,ridx,query,reference,tani,gani,ani,qcov,rcov,num_alns,len_ratio",
        "--threads",
        "0",
        "--verbose",
        "2",
    ]
    
    # Check for prefiltered pairs from kmer-db
    doc_dir = os.path.dirname(settings.doc_paths.matrix)
    filter_file = os.path.join(doc_dir, "kmerdb_filter.txt")
    if os.path.exists(filter_file):
        # Use --flt-kmerdb with a threshold of 0.0 to include all pairs in the file
        cmd.extend(["--flt-kmerdb", filter_file, "0.0"])
        print(f"Using kmer-db prefilter from: {filter_file}")
    # Not actual alignments, but can give insights into the analysis
    if settings.export_alignments and settings.alignment_export_path:
        alignment_file = os.path.join(
            settings.alignment_export_path, "lzani_alignments.txt"
        )
        cmd.extend(["--out-alignment", alignment_file])

        os.makedirs(settings.alignment_export_path, exist_ok=True)

    if settings.lzani.aw is not None:
        cmd.extend(["--aw", str(settings.lzani.aw)])

    if settings.lzani.am is not None:
        cmd.extend(["--am", str(settings.lzani.am)])

    if settings.lzani.mal is not None:
        cmd.extend(["--mal", str(settings.lzani.mal)])

    if settings.lzani.msl is not None:
        cmd.extend(["--msl", str(settings.lzani.msl)])

    if settings.lzani.mrd is not None:
        cmd.extend(["--mrd", str(settings.lzani.mrd)])

    if settings.lzani.mqd is not None:
        cmd.extend(["--mqd", str(settings.lzani.mqd)])

    if settings.lzani.reg is not None:
        cmd.extend(["--reg", str(settings.lzani.reg)])

    if settings.lzani.ar is not None:
        cmd.extend(["--ar", str(settings.lzani.ar)])

    process = subprocess.Popen(
        cmd,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,  # Merge stderr into stdout
        text=True,
        bufsize=1,  # Line buffered
    )

    output_lines = []
    pairs_pattern = re.compile(r"Pairs:\s*(\d+(?:\.\d+)?)%")

    try:
        if process.stdout is not None:
            for line in process.stdout:
                if canceled.value:
                    process.kill()
                    process.wait()
                    return result._replace(error="LZ-ANI run was canceled")

                output_lines.append(line)
                print(line, end="")

                # Extract pairs percentage
                match = pairs_pattern.search(line)
                if match:
                    progress = int(float(match.group(1)))
                    set_progress(progress)

        process.wait(timeout=300)
    except subprocess.TimeoutExpired:
        process.kill()
        process.wait()
        return result._replace(error="LZ-ANI timed out after 5 minutes")

    if process.returncode != 0:
        return result._replace(error=f"Error running LZ-ANI: {''.join(output_lines)}")


    # Check sequence count from IDs file
    import pandas as pd
    ids_df = pd.read_csv(settings.doc_paths.lzani_results_ids, sep="\t")
    sequence_count = len(ids_df)
    
    if sequence_count > 2500:
        # For large datasets, skip matrix building to save memory
        # UMAP will read TSV directly, other views are disabled
        print(f"Skipping matrix construction for large dataset ({sequence_count} sequences > 2500)")
        # Return minimal matrix to satisfy postprocessing
        import numpy as np
        ids = ids_df["id"].tolist()
        # Create a dummy 1x1 matrix to avoid errors in postprocessing
        dummy_matrix = np.zeros((1, 1))
        return result._replace(distance_matrix=dummy_matrix, ordered_ids=ids)
    else:
        # For smaller datasets, build the full matrix as usual
        matrix, ids = lzani_tsv_to_distance_matrix(
            settings.doc_paths.lzani_results,
            settings.doc_paths.lzani_results_ids,
            settings.lzani.score_type,
        )
        return result._replace(distance_matrix=matrix, ordered_ids=ids)
