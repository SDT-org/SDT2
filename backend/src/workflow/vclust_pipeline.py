import os
import logging
import tempfile
import subprocess
import shutil
from pathlib import Path
import sys
from typing import Dict, List, Tuple, Callable, Optional, Set
from multiprocessing.sharedctypes import Synchronized

import pandas as pd
from scipy.sparse import coo_matrix

from workflow.models import RunSettings, WorkflowResult, VclustSettings
from app_state import update_document

logger = logging.getLogger(__name__)

# Get the path to the original vclust script
# Find the root directory and then the vclust directory
current_file = Path(__file__).resolve()
root_dir = current_file.parent.parent.parent.parent
VCLUST_DIR = root_dir / "vclust"
VCLUST_SCRIPT = VCLUST_DIR / "vclust.py"

# Log the script path for debugging
logger.info(f"Using vclust script at: {VCLUST_SCRIPT}")

# Set paths to the binaries that vclust needs
KMERDB_BIN = root_dir / "3rd_party" / "kmer-db" / "bin" / "kmer-db.exe"
LZANI_BIN = root_dir / "3rd_party" / "lz-ani" / "bin" / "lz-ani.exe"
CLUSTY_BIN = root_dir / "3rd_party" / "clusty" / "bin" / "clusty.exe" 
MFASTA_BIN = root_dir / "3rd_party" / "ref-utils" / "bin" / "mfasta-tool.exe"

# Log the binary paths for debugging
logger.info(f"Using kmer-db at: {KMERDB_BIN}")
logger.info(f"Using lz-ani at: {LZANI_BIN}")
logger.info(f"Using clusty at: {CLUSTY_BIN}")
logger.info(f"Using mfasta-tool at: {MFASTA_BIN}")


def run(
    result: WorkflowResult,
    settings: RunSettings,
    set_progress: Callable[[int], None],
    canceled: Synchronized,
) -> WorkflowResult:
    """Run the vclust pipeline using the original vclust.py script."""
    try:
        # All intermediate files will be created in a temporary directory,
        # which is automatically cleaned up when the pipeline is finished.
        with tempfile.TemporaryDirectory() as temp_dir:
            temp_path = Path(temp_dir)
            logger.info(f"Created temporary directory for vclust intermediates: {temp_path}")
            
            # Check if vclust script exists
            if not VCLUST_SCRIPT.exists():
                return result._replace(error=f"vclust script not found at {VCLUST_SCRIPT}")
            
            # Setup paths for intermediate files
            prefilter_path = temp_path / "prefilter.txt"
            align_path = temp_path / "align.tsv"
            ids_path = temp_path / "ids.tsv"
            cluster_path = temp_path / "clusters.tsv"
            
            # Get vclust settings
            vclust_settings = settings.vclust or VclustSettings()
            
            # Step 1: Prefilter (30% progress)
            set_progress(0)
            logger.info("Step 1: Prefiltering sequences with kmer-db")
            
            prefilter_cmd = [
                sys.executable,
                str(VCLUST_SCRIPT),
                "prefilter",
                "-i", settings.fasta_path,
                "-o", str(prefilter_path),
                "--min-kmers", str(vclust_settings.kmer_min_kmers),
                "--min-ident", str(vclust_settings.kmer_min_similarity),
                "--kmers-fraction", str(vclust_settings.kmer_fraction),
                "-t", str(vclust_settings.threads if vclust_settings.threads > 0 else 0),
                "-v", str(vclust_settings.verbosity_level) if hasattr(vclust_settings, "verbosity_level") else "1"
            ]
            
            # Log the full command for debugging
            command_str = " ".join(prefilter_cmd)
            logger.info(f"Executing prefilter command: {command_str}")
            
            process = subprocess.Popen(
                prefilter_cmd,
                stdout=subprocess.PIPE,
                stderr=subprocess.STDOUT,  # Merge stderr with stdout
                text=True,
                bufsize=1
            )
            
            output_lines = []
            if process.stdout:
                for line in process.stdout:
                    line = line.strip()
                    output_lines.append(line)
                    logger.debug(line)
                    if "%" in line:
                        try:
                            progress_str = line.split("%")[0]
                            progress = int(float(progress_str))
                            set_progress(int(progress * 0.3)) # Prefilter is 30% of the work
                        except (ValueError, IndexError):
                            pass
                    if canceled.value:
                        process.kill()
                        return result._replace(error="Vclust pipeline was canceled")
            
            process.wait()
            if process.returncode != 0:
                error_msg = f"Prefilter step failed with return code {process.returncode}"
                error_details = "\n".join(output_lines)
                error_msg += f"\nError details: {error_details}"
                return result._replace(error=error_msg)
            
            set_progress(30)
            
            # Step 2: Align (60% progress)
            logger.info("Step 2: Aligning filtered sequence pairs with LZ-ANI")
            
            # Determine the score type (ani, tani, gani)
            score_type = getattr(settings.lzani, "score_type", "tani")
            
            align_cmd = [
                sys.executable,
                str(VCLUST_SCRIPT),
                "align",
                "-i", settings.fasta_path,
                "-o", str(align_path),
                "--filter", str(prefilter_path),
                "--outfmt", "complete",  # Get all fields for maximum information
                "--out-" + score_type, "0",  # No filtering at this stage
                "-t", str(vclust_settings.threads if vclust_settings.threads > 0 else 0),
                "-v", str(vclust_settings.verbosity_level) if hasattr(vclust_settings, "verbosity_level") else "1"
            ]
            
            # Add LZ-ANI parameters if provided in vclust_settings
            if hasattr(vclust_settings, "lzani_mal"):
                align_cmd.extend(["--mal", str(vclust_settings.lzani_mal)])
            if hasattr(vclust_settings, "lzani_msl"):
                align_cmd.extend(["--msl", str(vclust_settings.lzani_msl)])
            if hasattr(vclust_settings, "lzani_mrd"):
                align_cmd.extend(["--mrd", str(vclust_settings.lzani_mrd)])
            if hasattr(vclust_settings, "lzani_mqd"):
                align_cmd.extend(["--mqd", str(vclust_settings.lzani_mqd)])
            if hasattr(vclust_settings, "lzani_reg"):
                align_cmd.extend(["--reg", str(vclust_settings.lzani_reg)])
            if hasattr(vclust_settings, "lzani_aw"):
                align_cmd.extend(["--aw", str(vclust_settings.lzani_aw)])
            if hasattr(vclust_settings, "lzani_am"):
                align_cmd.extend(["--am", str(vclust_settings.lzani_am)])
            if hasattr(vclust_settings, "lzani_ar"):
                align_cmd.extend(["--ar", str(vclust_settings.lzani_ar)])
            
            # Log the full command for debugging
            command_str = " ".join(align_cmd)
            logger.info(f"Executing align command: {command_str}")
            
            process = subprocess.Popen(
                align_cmd,
                stdout=subprocess.PIPE,
                stderr=subprocess.STDOUT,  # Merge stderr with stdout
                text=True,
                bufsize=1
            )
            
            output_lines = []
            if process.stdout:
                for line in process.stdout:
                    line = line.strip()
                    output_lines.append(line)
                    logger.debug(line)
                    if "%" in line:
                        try:
                            # Extract progress percentage
                            progress_str = line.split("%")[0]
                            progress = int(float(progress_str))
                            set_progress(30 + int(progress * 0.6))
                        except (ValueError, IndexError):
                            pass
                    
                    if canceled.value:
                        process.kill()
                        return result._replace(error="Vclust pipeline was canceled")
            
            process.wait()
            if process.returncode != 0:
                error_msg = f"Align step failed with return code {process.returncode}"
                error_details = "\n".join(output_lines)
                error_msg += f"\nError details: {error_details}"
                return result._replace(error=error_msg)
            
            # Copy the output files to the document paths
            os.makedirs(os.path.dirname(settings.doc_paths.lzani_results), exist_ok=True)
            
            # Copy align.tsv to lzani_results
            if not os.path.exists(settings.doc_paths.lzani_results) or not os.path.samefile(align_path, settings.doc_paths.lzani_results):
                shutil.copy2(align_path, settings.doc_paths.lzani_results)
            
            # Get the IDs file
            # Note: vclust creates an IDs file alongside the align output but with .ids.tsv extension
            ids_file = align_path.with_suffix('.ids.tsv')
            if ids_file.exists():
                if not os.path.exists(settings.doc_paths.lzani_results_ids) or not os.path.samefile(ids_file, settings.doc_paths.lzani_results_ids):
                    shutil.copy2(ids_file, settings.doc_paths.lzani_results_ids)
            else:
                return result._replace(error="IDs file not found after alignment")
            
            set_progress(90)
            
            # Step 3: Cluster (10% progress)
            logger.info("Step 3: Clustering sequences")

            # Check if the align and ids files are valid before proceeding
            if not align_path.exists() or align_path.stat().st_size == 0:
                return result._replace(error="Align results file is empty or not found.")
            if not ids_file.exists() or ids_file.stat().st_size == 0:
                return result._replace(error="IDs file is empty or not found.")
            
            cluster_algorithm = getattr(vclust_settings, "cluster_algorithm", "cd-hit")
            cluster_metric = getattr(vclust_settings, "cluster_metric", "tani")
            cluster_threshold = getattr(vclust_settings, "cluster_threshold", 0.7)
            use_representatives = getattr(vclust_settings, "use_representatives", True)
            
            cluster_cmd = [
                sys.executable,
                str(VCLUST_SCRIPT),
                "cluster",
                "-i", settings.doc_paths.lzani_results,
                "--ids", settings.doc_paths.lzani_results_ids,
                "-o", str(cluster_path),
                "--algorithm", cluster_algorithm,
                "--metric", cluster_metric,
                f"--{cluster_metric}", str(cluster_threshold),
            ]
            
            # Add optional cluster parameters
            if use_representatives:
                cluster_cmd.append("-r")
            
            min_coverage = getattr(vclust_settings, "min_coverage", 0)
            if min_coverage > 0:
                cluster_cmd.extend(["--qcov", str(min_coverage), "--rcov", str(min_coverage)])
            
            min_length_ratio = getattr(vclust_settings, "min_length_ratio", 0)
            if min_length_ratio > 0:
                cluster_cmd.extend(["--len_ratio", str(min_length_ratio)])
            
            max_alignments = getattr(vclust_settings, "max_alignments", 0)
            if max_alignments > 0:
                cluster_cmd.extend(["--num_alns", str(max_alignments)])
            
            # Log the full command for debugging
            command_str = " ".join(cluster_cmd)
            logger.info(f"Executing cluster command: {command_str}")
            print(f"Executing cluster command: {command_str}") # Also print to stdout for visibility
            
            process = subprocess.Popen(
                cluster_cmd,
                stdout=subprocess.PIPE,
                stderr=subprocess.STDOUT,  # Merge stderr with stdout
                text=True,
                bufsize=1
            )
            
            output_lines = []
            if process.stdout:
                for line in process.stdout:
                    line = line.strip()
                    output_lines.append(line)
                    logger.debug(line)
                    if canceled.value:
                        process.kill()
                        return result._replace(error="Vclust pipeline was canceled")
            
            process.wait()
            if process.returncode != 0:
                error_msg = f"Cluster step failed with return code {process.returncode}"
                error_details = "\n".join(output_lines)
                error_msg += f"\nError details: {error_details}"
                return result._replace(error=error_msg)
            
            # Parse the clustering results to get representative sequences
            representatives = set()
            if cluster_path.exists():
                with open(cluster_path, 'r') as f:
                    for line in f:
                        parts = line.strip().split('\t')
                        if len(parts) >= 2:
                            seq_id = parts[0]
                            # If using representatives and format includes a representative flag
                            if use_representatives and len(parts) >= 3:
                                is_representative = parts[2] == "1"
                                if is_representative:
                                    representatives.add(seq_id)
                            else:
                                # Just add all sequences - they will be filtered later if needed
                                representatives.add(seq_id)
            
            # Create sparse distance matrix for UMAP
            sparse_matrix, sequence_ids = create_sparse_matrix(
                align_path=Path(settings.doc_paths.lzani_results),
                ids_path=Path(settings.doc_paths.lzani_results_ids),
                representative_ids=representatives,
                score_type=score_type
            )
            
            result = result._replace(
                distance_matrix=sparse_matrix,
                ordered_ids=sequence_ids,
            )
            
            # For large datasets, the only valid view is UMAP
            if len(sequence_ids) > 1000:
                doc_id = os.path.basename(os.path.dirname(settings.doc_paths.matrix))
                update_document(doc_id, dataView='umap')
            
            set_progress(100)
            logger.info(f"Pipeline complete: {len(representatives)} representative sequences found")
            
            return result
            
    except Exception as e:
        logger.error(f"Vclust pipeline error: {str(e)}", exc_info=True)
        return result._replace(error=f"Vclust pipeline failed: {str(e)}")


def create_sparse_matrix(
    align_path: Path,
    ids_path: Path,
    representative_ids: Set[str],
    score_type: str
) -> Tuple[coo_matrix, List[str]]:
    """Create a sparse distance matrix from alignment results for UMAP visualization."""
    # Read sequence IDs
    ids_df = pd.read_csv(ids_path, sep="\t")
    sequence_ids = ids_df["id"].tolist()
    id_to_idx = {seq_id: i for i, seq_id in enumerate(sequence_ids)}
    n_sequences = len(sequence_ids)
    
    row_indices, col_indices, distances = [], [], []
    
    # Read alignment results
    if align_path.exists() and align_path.stat().st_size > 0:
        try:
            # Read in chunks to handle large files
            for chunk in pd.read_csv(align_path, sep="\t", chunksize=10000):
                # Filter to include only representatives if specified
                if representative_ids:
                    chunk = chunk[
                        chunk['query'].isin(representative_ids) & 
                        chunk['reference'].isin(representative_ids)
                    ]
                
                if not chunk.empty:
                    # Calculate distance from similarity
                    chunk['distance'] = 1.0 - chunk[score_type]
                    
                    # Map sequence IDs to indices
                    chunk['i'] = chunk['query'].map(id_to_idx)
                    chunk['j'] = chunk['reference'].map(id_to_idx)
                    
                    # Filter out rows with missing indices
                    valid_mask = chunk['i'].notna() & chunk['j'].notna()
                    chunk = chunk[valid_mask]
                    
                    # Add to COO matrix data
                    row_indices.extend(chunk['i'].astype(int).tolist())
                    col_indices.extend(chunk['j'].astype(int).tolist())
                    distances.extend(chunk['distance'].tolist())
                    
                    # Add symmetric entries (i,j) -> (j,i)
                    row_indices.extend(chunk['j'].astype(int).tolist())
                    col_indices.extend(chunk['i'].astype(int).tolist())
                    distances.extend(chunk['distance'].tolist())
        except Exception as e:
            logger.error(f"Error creating sparse matrix: {str(e)}", exc_info=True)
    
    # If no valid pairs, create empty matrix
    if not distances:
        logger.warning("No valid pairs found in alignment results. UMAP may fail.")
        sparse_matrix = coo_matrix((n_sequences, n_sequences))
    else:
        sparse_matrix = coo_matrix(
            (distances, (row_indices, col_indices)),
            shape=(n_sequences, n_sequences)
        )
    
    return sparse_matrix, sequence_ids


def estimate_vclust_memory(
    num_sequences: int,
    avg_sequence_length: int,
    kmer_fraction: float = 0.2
) -> Dict[str, float]:
    """Estimate memory usage for the vclust pipeline."""
    # Kmer-db memory estimate
    base_memory_per_seq = 0.001  # GB
    kmer_memory_per_seq = (avg_sequence_length / 1000) * 0.0001 * kmer_fraction
    kmerdb_memory = num_sequences * (base_memory_per_seq + kmer_memory_per_seq)
    
    # Apply overhead factor for sparse matrix
    kmerdb_memory *= 1.2 if num_sequences > 10000 else 1.5
    
    # Estimate number of pairs (assume 2% of all possible pairs pass prefilter)
    estimated_pairs = int(num_sequences * num_sequences * 0.02)
    
    # LZ-ANI memory estimate (per pair)
    lzani_memory = estimated_pairs * 0.000001  # GB
    
    # Clusty memory estimate
    clusty_memory = estimated_pairs * 0.000002  # GB
    
    return {
        "kmerdb": kmerdb_memory,
        "lzani": lzani_memory,
        "clusty": clusty_memory,
        "total": kmerdb_memory + lzani_memory + clusty_memory
    }
