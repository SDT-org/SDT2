import os
import subprocess
import tempfile
from pathlib import Path
from typing import List, Tuple, Optional, Dict, Set
import logging

logger = logging.getLogger(__name__)


def run_clusty(
    ani_file: Path,
    ids_file: Path,
    output_path: Path,
    algorithm: str = "single",
    threshold: float = 0.95,
    metric: str = "tani",
    min_ani: Optional[float] = None,
    min_af: Optional[float] = None,
    threads: Optional[int] = None
) -> Dict[str, List[str]]:
    if not ani_file.exists():
        raise FileNotFoundError(f"ANI file not found: {ani_file}")
    if not ids_file.exists():
        raise FileNotFoundError(f"IDs file not found: {ids_file}")
    
    clusty_exe = Path(__file__).parent.parent.parent.parent / "bin" / "clusty.exe"
    if not clusty_exe.exists():
        raise FileNotFoundError(f"clusty executable not found at {clusty_exe}")
    
    if threads is None:
        threads = os.cpu_count() or 1
    
    output_path.parent.mkdir(parents=True, exist_ok=True)

    # Build the command using the correct argument format for clusty.exe
    cmd = [
        str(clusty_exe),
        '--objects-file', str(ids_file),
        '--algo', algorithm,
        '--id-cols', 'qidx', 'ridx',
        '--distance-col', metric,
        '--similarity',
        '--numeric-ids',
        str(ani_file),
        str(output_path)
    ]

    # Add filters. Clusty expects the threshold to be passed as a --min filter.
    cmd.extend(['--min', metric, str(threshold)])

    if min_ani is not None:
        cmd.extend(['--min', 'ani', str(min_ani)])
    if min_af is not None:
        # Assuming min_af corresponds to coverage. Clusty uses qcov and rcov.
        # We will apply it to both for simplicity.
        cmd.extend(['--min', 'qcov', str(min_af)])
        cmd.extend(['--min', 'rcov', str(min_af)])
    
    logger.info(f"Running Clusty {algorithm} clustering with threshold={threshold}")
    result = subprocess.run(cmd, capture_output=True, text=True)
    
    if result.returncode != 0:
        raise RuntimeError(f"Clusty failed: {result.stderr}")
    
    # Parse clustering results
    clusters = {}
    with open(output_path, 'r') as f:
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) >= 2:
                seq_id = parts[0]
                cluster_id = parts[1]
                
                if cluster_id not in clusters:
                    clusters[cluster_id] = []
                clusters[cluster_id].append(seq_id)
    
    logger.info(f"Clustering complete: {len(clusters)} clusters formed")
    return clusters


def get_cluster_representatives(
    clusters: Dict[str, List[str]],
    sequence_lengths: Dict[str, int]
) -> Dict[str, str]:
    representatives = {}
    
    for cluster_id, members in clusters.items():
        # Find longest sequence as representative
        longest_seq = max(members, key=lambda x: sequence_lengths.get(x, 0))
        representatives[cluster_id] = longest_seq
    
    return representatives


def calculate_votu_sizes(clusters: Dict[str, List[str]]) -> Dict[str, int]:
    return {cluster_id: len(members) for cluster_id, members in clusters.items()}


def filter_clusters_by_size(
    clusters: Dict[str, List[str]],
    min_size: int = 1,
    max_size: Optional[int] = None
) -> Dict[str, List[str]]:
    filtered = {}
    
    for cluster_id, members in clusters.items():
        size = len(members)
        if size >= min_size and (max_size is None or size <= max_size):
            filtered[cluster_id] = members
    
    return filtered


def filter_ani_results(
    ani_file: Path,
    ids_to_keep: Set[str],
    output_path: Path
):
    """
    Filters an ANI results file to include only pairs where both query and target
    are in the provided set of IDs.
    """
    logger.info(f"Filtering ANI results in {ani_file} to {len(ids_to_keep)} sequences.")
    with open(ani_file, 'r') as infile, open(output_path, 'w') as outfile:
        header = infile.readline()
        outfile.write(header)
        
        for line in infile:
            parts = line.strip().split('\t')
            if len(parts) >= 2:
                query_id, target_id = parts[0], parts[1]
                if query_id in ids_to_keep and target_id in ids_to_keep:
                    outfile.write(line)
    logger.info(f"Filtered ANI results saved to {output_path}")


def deduplicate_sequences(
    ani_file: Path,
    ids_file: Path,
    output_path: Path,
    threshold: float = 0.99,
    threads: Optional[int] = None
) -> Set[str]:
    clusters = run_clusty(
        ani_file=ani_file,
        ids_file=ids_file,
        output_path=output_path,
        algorithm="cd-hit",
        threshold=threshold,
        metric="tani",
        threads=threads
    )
    
    # Get sequence lengths from IDs file
    seq_lengths = {}
    with open(ids_file, 'r') as f:
        header = f.readline()
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) >= 3:
                seq_id = parts[0]
                length = int(parts[2])
                seq_lengths[seq_id] = length
    
    # If clustering is empty, it means no sequences met the cd-hit threshold.
    # In this case, we should return all the original sequences to allow the
    # pipeline to continue, rather than causing it to fail.
    if not clusters:
        logger.warning(
            f"No sequences met the cd-hit threshold of {threshold}. "
            f"Skipping deduplication and passing all sequences to the next step."
        )
        return set(seq_lengths.keys())

    representatives = get_cluster_representatives(clusters, seq_lengths)
    return set(representatives.values())
