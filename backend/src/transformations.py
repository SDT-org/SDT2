from numpy import tril, triu_indices, around, nan
import numpy as np
import pandas as pd
from pandas import DataFrame
import multiprocessing as mp
import threading
from concurrent.futures import ThreadPoolExecutor
import os
import psutil


def to_triangle(matrix, convert_to_similarity=True, fill_value=nan or None):

    is_dataframe = isinstance(matrix, DataFrame)

    # Get numpy array and metadata
    if is_dataframe:
        data = matrix.to_numpy()
        index = matrix.index
        columns = matrix.columns
    else:
        data = matrix

    # Convert to similarity
    if convert_to_similarity:
        data = 100 - data

    # Create lower triangle
    result = tril(around(data, 2))

    # Fill upper triangle
    upper_indices = triu_indices(result.shape[0], k=1)
    if fill_value is None:
        result[upper_indices] = None
        result = np.where(np.isnan(result), None, result)  # type: ignore[arg-type]
    else:
        result[upper_indices] = fill_value

    # Return in original format
    if is_dataframe:
        return DataFrame(result, index=index, columns=columns)  # type: ignore[arg-type]
    else:
        return result


def similarity_triangle_to_matrix(matrix_lower: DataFrame) -> DataFrame:
    index = matrix_lower.index
    sim_lower = np.round(np.tril(100 - matrix_lower), 2)
    df = sim_lower + np.triu(sim_lower.T, 1)
    return DataFrame(df, index=index)


def read_csv_file(
    filepath, sep=",", header=None, index_col=None, as_list=False, extract_column=None
):
    # Special handling for matrix files with variable columns
    if header is None and index_col == 0 and sep == ",":
        with open(filepath, "r") as temp_f:
            col_count = [len(l.split(sep)) for l in temp_f.readlines()]
            column_names = [i for i in range(0, max(col_count))]
        df = pd.read_csv(
            filepath,
            delimiter=sep,
            index_col=index_col,
            header=header,
            names=column_names,
        )
    else:
        df = pd.read_csv(filepath, sep=sep, header=header, index_col=index_col)

    # Return based on options
    if extract_column:
        return df[extract_column].tolist()
    elif as_list:
        return df.values.tolist()
    else:
        return df


def read_csv_matrix(filepath):
    return read_csv_file(filepath, sep=",", header=None, index_col=0)


def read_stats_csv(filepath):
    return read_csv_file(filepath, header=0)


def read_columns_csv(filepath):
    # need typimg for id coliumn  **DtypeWarning: Columns (2) have mixed types. Specify dtype option on import or set low_memory=False.
    dtype_spec = {
        0: str,
        1: str,
        2: float,
    }  # First Sequence, Second Sequence, Identity Score
    df = pd.read_csv(filepath, header=0, dtype=dtype_spec)
    return df.values.tolist()


def read_tsv_file(filepath, id_column=None):
    return read_csv_file(filepath, sep="\t", extract_column=id_column)


def _process_chunk_vectorized_memory_efficient(chunk, id_to_idx, score_column):
    """Memory-efficient vectorized chunk processing that returns sparse data."""
    if chunk.empty:
        return [], [], []

    # Vectorized operations instead of iterrows()
    query_indices = chunk["query"].map(id_to_idx).fillna(-1).astype(int)
    ref_indices = chunk["reference"].map(id_to_idx).fillna(-1).astype(int)

    # Filter out invalid mappings
    valid_mask = (query_indices >= 0) & (ref_indices >= 0)
    query_valid = query_indices[valid_mask].values
    ref_valid = ref_indices[valid_mask].values
    scores_valid = (chunk[score_column][valid_mask] * 100).values

    return query_valid, ref_valid, scores_valid


def lzani_tsv_to_distance_matrix(
    results_tsv_path, ids_tsv_path, score_column="ani", chunksize=10000, n_workers=None
):
    print(f"Building distance matrix from {os.path.basename(results_tsv_path)}...")

    # Read IDs - explicitly get list
    ids_df = pd.read_csv(ids_tsv_path, sep="\t")
    all_ids = ids_df["id"].tolist()
    n_ids = len(all_ids)
    print(f"Processing {n_ids} sequences")

    # Show memory usage
    if psutil:
        process = psutil.Process()
        memory_mb = process.memory_info().rss / 1024 / 1024
        print(f"Current memory usage: {memory_mb:.1f} MB")

    # Create ID to index mapping for fast lookup
    id_to_idx = {id_val: i for i, id_val in enumerate(all_ids)}

    # Initialize shared matrix directly
    matrix_size_mb = (n_ids * n_ids * 8) / 1024 / 1024  # 8 bytes per float64
    print(f"Allocating {n_ids}x{n_ids} matrix ({matrix_size_mb:.1f} MB)")
    matrix_np = np.zeros((n_ids, n_ids))

    # Determine number of workers
    if n_workers is None:
        n_workers = min(mp.cpu_count(), 8)  # Cap at 8 to avoid too many threads

    # Threading approach for memory efficiency
    try:
        print(f"Reading TSV file in chunks of {chunksize} rows...")
        # Pre-read all chunks to avoid concurrent file I/O issues
        chunk_reader = pd.read_csv(results_tsv_path, sep="\t", chunksize=chunksize)
        chunks = list(chunk_reader)

        if not chunks:
            print("No chunks found in TSV file")
        else:
            print(
                f"Loaded {len(chunks)} chunks, processing with {n_workers} threads..."
            )

            # Thread lock for safe matrix updates
            matrix_lock = threading.Lock()

            def process_chunk_worker(chunk):
                """Worker function that processes a chunk and updates shared matrix."""
                # Process chunk and get sparse representation
                query_indices, ref_indices, scores = (
                    _process_chunk_vectorized_memory_efficient(
                        chunk, id_to_idx, score_column
                    )
                )

                # Thread-safe matrix update
                if len(query_indices) > 0:
                    with matrix_lock:
                        matrix_np[query_indices, ref_indices] = scores

                return len(query_indices)  # Return number of entries processed

            # Process chunks in parallel using threading
            if n_workers > 1 and len(chunks) > 1:
                print("Starting threaded chunk processing...")
                with ThreadPoolExecutor(max_workers=n_workers) as executor:
                    # Submit all chunks and track progress
                    futures = [
                        executor.submit(process_chunk_worker, chunk) for chunk in chunks
                    ]

                    # Wait for completion and track progress
                    completed = 0
                    for future in futures:
                        future.result()
                        completed += 1
                        if completed % 10 == 0:
                            print(f"Completed {completed}/{len(chunks)} chunks")

                print(f"Threaded processing completed: {len(chunks)} chunks")
            else:
                # Sequential fallback for small datasets
                print("Using sequential processing...")
                for i, chunk in enumerate(chunks):
                    process_chunk_worker(chunk)
                    if (i + 1) % 10 == 0:
                        print(f"Processed chunk {i+1}/{len(chunks)}")

    except pd.errors.EmptyDataError:
        print("Warning: TSV file appears to be empty")
        pass

    # Replace exact 0s (unaligned) with 0 - no need for NaN handling since we init with zeros
    matrix_np = np.where(matrix_np == 0, 0, matrix_np)

    # LZANI has slightly different scores for query v. reference vs. reference v query
    # Only average if both directions have valid alignments (>1% threshold)
    print("Computing matrix symmetry (averaging bidirectional scores)...")

    # Memory-efficient in-place symmetry calculation
    print(
        f"Using memory-efficient in-place symmetry calculation for {n_ids}x{n_ids} matrix..."
    )

    # Process in chunks to avoid memory explosion
    chunk_size = min(500, n_ids)  # Process 500 rows at a time

    for start_i in range(0, n_ids, chunk_size):
        end_i = min(start_i + chunk_size, n_ids)
        if start_i % (chunk_size * 10) == 0:
            print(f"Processing symmetry for rows {start_i}-{end_i-1}/{n_ids}")

        for i in range(start_i, end_i):
            # Set diagonal to 100 (self-comparison)
            matrix_np[i, i] = 100.0

            # Process upper triangle for this row
            for j in range(i + 1, n_ids):
                forward = matrix_np[i, j]
                reverse = matrix_np[j, i]

                # Calculate symmetric value
                if forward > 1 and reverse > 1:
                    sym_value = (forward + reverse) / 2
                elif forward > 1:
                    sym_value = forward
                elif reverse > 1:
                    sym_value = reverse
                else:
                    sym_value = 0

                # Set both symmetric positions
                matrix_np[i, j] = sym_value
                matrix_np[j, i] = sym_value

    print("Matrix symmetry calculation completed")

    # Convert similarity to distance
    print("Converting similarity scores to distance matrix...")
    distance_matrix = 100 - matrix_np

    print(f"Distance matrix construction completed ({n_ids}x{n_ids})")
    return distance_matrix, all_ids
