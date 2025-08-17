import os
import logging
import tempfile
from pathlib import Path
from typing import Optional, Dict, List, Tuple, Callable
from multiprocessing.sharedctypes import Synchronized

from workflow.models import RunSettings, WorkflowResult, VclustSettings
from workflow.analyze.kmerdb import run_kmerdb_prefilter
from workflow.cluster.clusty import deduplicate_sequences, filter_ani_results
import pandas as pd
from scipy.sparse import coo_matrix
from app_state import update_document

logger = logging.getLogger(__name__)


def run(
    result: WorkflowResult,
    settings: RunSettings,
    set_progress: Callable[[int], None],
    canceled: Synchronized,
) -> WorkflowResult:
    try:
        # All intermediate files will be created in a temporary directory,
        # which is automatically cleaned up when the pipeline is finished.
        with tempfile.TemporaryDirectory() as temp_dir:
            temp_path = Path(temp_dir)
            logger.info(f"Created temporary directory for pipeline intermediates: {temp_path}")

            # Stage 1: K-mer prefiltering (30% of progress)
            set_progress(0)
            logger.info("Stage 1: K-mer prefiltering")
            
            filter_path = temp_path / "kmerdb_filter.txt"
            
            # Use configurable settings or defaults
            vclust_settings = settings.vclust or VclustSettings()
            kmer_result = run_kmerdb_prefilter(
                fasta_path=Path(settings.fasta_path),
                output_path=filter_path,
                min_similarity=vclust_settings.kmer_min_similarity,
                min_kmers=vclust_settings.kmer_min_kmers,
                kmer_fraction=vclust_settings.kmer_fraction,
                threads=1
            )
            
            if canceled.value:
                return result._replace(error="Vclust pipeline was canceled")
            
            set_progress(30)
            
            # Stage 2: Run LZANI on filtered pairs (60% of progress)
            logger.info("Stage 2: Running LZ-ANI on filtered pairs")
            
            from workflow.analyze import lzani
            result = lzani.run(result, settings, lambda p: set_progress(30 + int(p * 0.6)), canceled)
            
            if result.error:
                return result
            
            set_progress(90)
            
            # Stage 3: Deduplication with cd-hit (10% of progress)
            logger.info("Stage 3: Deduplicating sequences with cd-hit")
            
            dedup_cluster_path = temp_path / "cdhit_clusters.txt"
            representative_ids = deduplicate_sequences(
                ani_file=Path(settings.doc_paths.lzani_results),
                ids_file=Path(settings.doc_paths.lzani_results_ids),
                output_path=dedup_cluster_path,
                threshold=vclust_settings.cdhit_threshold,
                threads=None
            )
            
            logger.info(f"Reduced to {len(representative_ids)} representative sequences.")
            
            # Filter the ANI results to only include representatives
            filtered_ani_path = temp_path / "lzani_results_filtered.tsv"
            filter_ani_results(
                ani_file=Path(settings.doc_paths.lzani_results),
                ids_to_keep=representative_ids,
                output_path=filtered_ani_path
            )
            
            # After filtering, create a sparse distance matrix for the UMAP step.
            logger.info("Creating sparse distance matrix from filtered ANI results for UMAP.")
            
            ids_df = pd.read_csv(settings.doc_paths.lzani_results_ids, sep="\t")
            sequence_ids = ids_df["id"].tolist()
            id_to_idx = {seq_id: i for i, seq_id in enumerate(sequence_ids)}
            n_sequences = len(sequence_ids)

            row_indices, col_indices, distances = [], [], []
            score_type = getattr(settings.lzani, 'score_type', 'tani')

            if filtered_ani_path.exists() and filtered_ani_path.stat().st_size > 0:
                ani_df = pd.read_csv(filtered_ani_path, sep="\t")
                if not ani_df.empty:
                    ani_df['distance'] = 1.0 - ani_df[score_type]
                    ani_df['i'] = ani_df['query'].map(id_to_idx)
                    ani_df['j'] = ani_df['reference'].map(id_to_idx)
                    
                    valid_mask = ani_df['i'].notna() & ani_df['j'].notna()
                    ani_df = ani_df[valid_mask]

                    row_indices.extend(ani_df['i'].astype(int).tolist())
                    col_indices.extend(ani_df['j'].astype(int).tolist())
                    distances.extend(ani_df['distance'].tolist())
                    
                    row_indices.extend(ani_df['j'].astype(int).tolist())
                    col_indices.extend(ani_df['i'].astype(int).tolist())
                    distances.extend(ani_df['distance'].tolist())

            if not distances:
                 logger.warning("No valid pairs found in filtered ANI results. UMAP may fail.")
                 sparse_matrix = coo_matrix((n_sequences, n_sequences))
            else:
                sparse_matrix = coo_matrix(
                    (distances, (row_indices, col_indices)),
                    shape=(n_sequences, n_sequences)
                )

            result = result._replace(
                distance_matrix=sparse_matrix,
                ordered_ids=sequence_ids,
            )
            
            # For large datasets, the only valid view is UMAP.
            # We must update the document state to switch the view automatically.
            if n_sequences > 1000:
                # Extract doc_id from the document path (it's the parent directory name)
                doc_id = os.path.basename(os.path.dirname(settings.doc_paths.matrix))
                update_document(doc_id, dataView='umap')

            set_progress(100)
            logger.info(f"Pipeline complete: Deduplication finished.")
            
            return result
            
    except Exception as e:
        logger.error(f"Vclust pipeline error: {str(e)}")
        return result._replace(error=f"Vclust pipeline failed: {str(e)}")


def estimate_vclust_memory(
    num_sequences: int,
    avg_sequence_length: int,
    kmer_fraction: float = 0.2
) -> Dict[str, float]:
    from workflow.analyze.kmerdb import estimate_memory_usage
    
    kmerdb_memory = estimate_memory_usage(num_sequences, avg_sequence_length, kmer_fraction)
    
    estimated_pairs = int(num_sequences * num_sequences * 0.02)
    lzani_memory = estimated_pairs * 0.000001
    
    clusty_memory = estimated_pairs * 0.000002
    
    return {
        "kmerdb": kmerdb_memory,
        "lzani": lzani_memory,
        "clusty": clusty_memory,
        "total": kmerdb_memory + lzani_memory + clusty_memory
    }
