import numpy as np
import pandas as pd
from hdbscan import HDBSCAN
from typing import Dict, List, Optional, Tuple


def run_hdbscan_clustering(
    distance_matrix: np.ndarray,
    sequence_ids: List[str],
    min_cluster_size: int = 5,
    min_samples: Optional[int] = None,
    cluster_selection_epsilon: float = 0.0,
    alpha: float = 1.0,
    metric: str = 'precomputed'
) -> Dict[str, int]:
    
    if min_samples is None:
        min_samples = min_cluster_size
    
    # Ensure the matrix is symmetric
    distance_matrix = (distance_matrix + distance_matrix.T) / 2
    
    # Run HDBSCAN
    clusterer = HDBSCAN(
        min_cluster_size=min_cluster_size,
        min_samples=min_samples,
        cluster_selection_epsilon=cluster_selection_epsilon,
        alpha=alpha,
        metric=metric,
        cluster_selection_method='eom'  # Excess of Mass
    )
    
    cluster_labels = clusterer.fit_predict(distance_matrix)
    
    # Create mapping from sequence ID to cluster
    # HDBSCAN uses -1 for noise points, we'll map these to cluster 0
    cluster_map = {}
    for i, seq_id in enumerate(sequence_ids):
        cluster_id = int(cluster_labels[i])  # Convert numpy int64 to Python int
        if cluster_id == -1:
            cluster_id = 0  # Noise/unclustered points
        else:
            cluster_id += 1  # Shift cluster IDs to start from 1
        cluster_map[seq_id] = int(cluster_id)  # Ensure it's a Python int
    
    return cluster_map


def get_hdbscan_clusters_dataframe(
    distance_matrix: np.ndarray,
    sequence_ids: List[str],
    min_cluster_size: int = 5,
    cluster_selection_epsilon: float = 0.0
) -> pd.DataFrame:
    """
    Get HDBSCAN clusters as a DataFrame.
    
    Returns:
    --------
    pd.DataFrame
        DataFrame with columns 'ID' and 'Cluster'
    """
    cluster_map = run_hdbscan_clustering(
        distance_matrix,
        sequence_ids,
        min_cluster_size=min_cluster_size,
        cluster_selection_epsilon=cluster_selection_epsilon
    )
    
    df = pd.DataFrame([
        {'ID': seq_id, 'Cluster': cluster_id}
        for seq_id, cluster_id in cluster_map.items()
    ])
    
    return df


def get_cluster_stats(cluster_map: Dict[str, int]) -> Dict:
    """
    Calculate statistics about the clustering.
    """
    cluster_sizes = {}
    for cluster_id in cluster_map.values():
        cluster_sizes[cluster_id] = cluster_sizes.get(cluster_id, 0) + 1
    
    # Separate noise (cluster 0) from actual clusters
    noise_count = cluster_sizes.get(0, 0)
    actual_clusters = {k: v for k, v in cluster_sizes.items() if k != 0}
    
    stats = {
        'total_sequences': int(len(cluster_map)),
        'total_clusters': int(len(actual_clusters)),
        'noise_points': int(noise_count),
        'largest_cluster_size': int(max(actual_clusters.values())) if actual_clusters else 0,
        'smallest_cluster_size': int(min(actual_clusters.values())) if actual_clusters else 0,
        'average_cluster_size': float(sum(actual_clusters.values()) / len(actual_clusters)) if actual_clusters else 0.0,
        'cluster_sizes': {int(k): int(v) for k, v in cluster_sizes.items()}
    }
    
    return stats
