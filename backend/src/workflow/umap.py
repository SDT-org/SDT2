import numpy as np
import pandas as pd
from umap import UMAP
import json
import os
import warnings
from typing import Dict, List, Optional, Tuple, Union
from dataclasses import dataclass, asdict
from transformations import read_csv_matrix
from scipy.sparse import coo_matrix, csr_matrix, issparse

# Suppress UMAP warnings
warnings.filterwarnings("ignore", category=UserWarning, module="umap")

# issues with numpy integers, had ai make oop classes for
@dataclass
class UMAPConfig:
    n_components: int = 2
    n_neighbors: int = 15
    min_dist: float = 0.1
    metric: str = 'precomputed'
    n_epochs: Optional[int] = None
    learning_rate: float = 1.0
    init: str = 'spectral'
    random_state: int = 42
    spread: float = 1.0
    negative_sample_rate: int = 5
    local_connectivity: float = 1.0
    repulsion_strength: float = 1.0
    
    def to_dict(self):
        return asdict(self)


@dataclass 
class UMAPResult:
    embedding: np.ndarray
    sequence_ids: List[str]
    config: UMAPConfig
    metadata: Optional[Dict] = None
    cluster_assignments: Optional[Dict[str, int]] = None
    
    def to_dataframe(self) -> pd.DataFrame:
        df = pd.DataFrame(
            self.embedding,
            columns=[f'UMAP{i+1}' for i in range(self.embedding.shape[1])]
        )
        df['sequence_id'] = self.sequence_ids
        
        if self.cluster_assignments:
            df['cluster'] = [self.cluster_assignments.get(sid, -1) for sid in self.sequence_ids]
            
        if self.metadata:
            for key, values in self.metadata.items():
                if isinstance(values, dict):
                    df[key] = [values.get(sid, '') for sid in self.sequence_ids]
                    
        return df


def run_umap(distance_matrix: Union[np.ndarray, coo_matrix, csr_matrix],
             sequence_ids: List[str],
             config: Optional[UMAPConfig] = None,
             cluster_assignments: Optional[Dict[str, int]] = None,
             metadata: Optional[Dict] = None,
             supervised: bool = False) -> UMAPResult:
    
    if distance_matrix.shape[0] != distance_matrix.shape[1]:
        raise ValueError("Distance matrix must be square")
    if len(sequence_ids) != distance_matrix.shape[0]:
        raise ValueError("Number of sequence IDs must match matrix dimensions")
    
    config = config or UMAPConfig()
    
    # Only symmetrize if it's a dense matrix
    if isinstance(distance_matrix, np.ndarray):
        distance_matrix = (distance_matrix + distance_matrix.T) / 2
    
    y = None
    if supervised and cluster_assignments:
        y = np.array([cluster_assignments.get(sid, -1) for sid in sequence_ids])
        valid_mask = y != -1
        if not all(valid_mask):
            distance_matrix = distance_matrix[valid_mask][:, valid_mask]
            sequence_ids = [sid for sid, valid in zip(sequence_ids, valid_mask) if valid]
            y = y[valid_mask]
    
    reducer = UMAP(**config.to_dict())
    embedding = reducer.fit_transform(distance_matrix, y) if y is not None else reducer.fit_transform(distance_matrix)

    
    return UMAPResult(
        embedding=embedding,
        sequence_ids=sequence_ids,
        config=config,
        metadata=metadata,
        cluster_assignments=cluster_assignments,

    )


def run_parameter_sweep(distance_matrix: np.ndarray,
                       sequence_ids: List[str],
                       n_neighbors_list: List[int] = [5, 15, 30, 50],
                       min_dist_list: List[float] = [0.0, 0.1, 0.25, 0.5],
                       cluster_assignments: Optional[Dict[str, int]] = None) -> Dict[Tuple[int, float], UMAPResult]:
    
    results = {}
    for n_neighbors in n_neighbors_list:
        for min_dist in min_dist_list:
            print(f"Running UMAP: n_neighbors={n_neighbors}, min_dist={min_dist}")
            config = UMAPConfig(n_neighbors=n_neighbors, min_dist=min_dist)
            results[(n_neighbors, min_dist)] = run_umap(
                distance_matrix, 
                sequence_ids,
                config=config,
                cluster_assignments=cluster_assignments
            )
    return results


def export_for_d3(result: UMAPResult, output_path: str):
    df = result.to_dataframe()
    
    data = {
        'embedding': [
            {
                'id': row['sequence_id'],
                'x': row['UMAP1'],
                'y': row['UMAP2'],
                'cluster': row.get('cluster', None),
                **{k: row.get(k, '') for k in result.metadata.keys() if result.metadata}
            }
            for _, row in df.iterrows()
        ],
        'config': result.config.to_dict(),
        'bounds': {
            'x': [float(df['UMAP1'].min()), float(df['UMAP1'].max())],
            'y': [float(df['UMAP2'].min()), float(df['UMAP2'].max())]
        }
    }
    
    with open(output_path, 'w') as f:
        json.dump(data, f, indent=2)


def batch_analysis(distance_matrix_path: str,
                  output_dir: str,
                  cluster_csv_paths: Dict[str, str] = None,
                  metadata_csv_path: str = None,
                  n_neighbors_list: List[int] = [5, 15, 30],
                  min_dist_list: List[float] = [0.1, 0.25, 0.5]):
    
    os.makedirs(output_dir, exist_ok=True)
    
    df_matrix = read_csv_matrix(distance_matrix_path)
    distance_matrix = df_matrix.to_numpy()
    sequence_ids = df_matrix.index.tolist()
    
    # Load metadata if provided
    metadata = None
    if metadata_csv_path and os.path.exists(metadata_csv_path):
        df_meta = pd.read_csv(metadata_csv_path)
        id_col = df_meta.columns[0]
        metadata = {
            col: dict(zip(df_meta[id_col], df_meta[col].astype(str)))
            for col in df_meta.columns if col != id_col
        }
    
    results = []
    
    # Run for each clustering method
    if cluster_csv_paths:
        for method_name, cluster_path in cluster_csv_paths.items():
            if os.path.exists(cluster_path):
                df_clusters = pd.read_csv(cluster_path)
                cluster_col = [col for col in df_clusters.columns if 'Cluster' in col][0]
                cluster_assignments = dict(zip(df_clusters['ID'], df_clusters[cluster_col]))
                
                # Run parameter sweep
                sweep_results = run_parameter_sweep(
                    distance_matrix,
                    sequence_ids,
                    n_neighbors_list,
                    min_dist_list,
                    cluster_assignments
                )
                
                for (n_neighbors, min_dist), result in sweep_results.items():
                    result.metadata = metadata
                    filename = f'umap_{method_name}_n{n_neighbors}_d{min_dist}.json'
                    export_for_d3(result, os.path.join(output_dir, filename))
                    results.append({
                        'method': method_name,
                        'n_neighbors': n_neighbors,
                        'min_dist': min_dist,
                        'file': filename
                    })
    else:
        # Run without clustering
        sweep_results = run_parameter_sweep(
            distance_matrix,
            sequence_ids,
            n_neighbors_list,
            min_dist_list
        )
        
        for (n_neighbors, min_dist), result in sweep_results.items():
            result.metadata = metadata
            filename = f'umap_n{n_neighbors}_d{min_dist}.json'
            export_for_d3(result, os.path.join(output_dir, filename))
            results.append({
                'n_neighbors': n_neighbors,
                'min_dist': min_dist,
                'file': filename
            })
    
    # Save manifest
    manifest = {
        'results': results,
        'distance_matrix': distance_matrix_path,
        'metadata': metadata_csv_path,
        'parameters': {
            'n_neighbors': n_neighbors_list,
            'min_dist': min_dist_list
        }
    }
    
    with open(os.path.join(output_dir, 'manifest.json'), 'w') as f:
        json.dump(manifest, f, indent=2)
    
    return results


def compare_clustering_methods(distance_matrix_path: str,
                             cluster_results_dir: str,
                             output_path: str,
                             n_neighbors: int = 15,
                             min_dist: float = 0.1):
    
    df_matrix = read_csv_matrix(distance_matrix_path)
    distance_matrix = df_matrix.to_numpy()
    sequence_ids = df_matrix.index.tolist()
    
    # Run UMAP once
    config = UMAPConfig(n_neighbors=n_neighbors, min_dist=min_dist)
    base_result = run_umap(distance_matrix, sequence_ids, config=config)
    
    comparison_data = {
        'embedding': base_result.embedding.tolist(),
        'sequence_ids': sequence_ids,
        'methods': {}
    }
    
    # Load clustering results from different methods
    for filename in os.listdir(cluster_results_dir):
        if filename.endswith('.csv') and 'cluster' in filename.lower():
            method_name = filename.replace('cluster_', '').replace('.csv', '')
            df_clusters = pd.read_csv(os.path.join(cluster_results_dir, filename))
            cluster_col = [col for col in df_clusters.columns if 'Cluster' in col][0]
            
            clusters = {}
            for _, row in df_clusters.iterrows():
                clusters[row['ID']] = int(row[cluster_col])
            
            comparison_data['methods'][method_name] = [
                clusters.get(sid, -1) for sid in sequence_ids
            ]
    
    with open(output_path, 'w') as f:
        json.dump(comparison_data, f, indent=2)
    
    return comparison_data


def run_umap_from_lzani_tsv(
    results_tsv_path: str,
    ids_tsv_path: str,
    config: Optional[UMAPConfig] = None,
    score_type: str = "ani",
    chunk_size: int = 100000  # Process TSV in chunks for large files
) -> UMAPResult:
    config = config or UMAPConfig()
    
    ids_df = pd.read_csv(ids_tsv_path, sep="\t")
    sequence_ids = ids_df["id"].tolist()
    n_sequences = len(sequence_ids)
    
    print(f"Processing UMAP for {n_sequences} sequences")
    
    id_to_idx = {seq_id: idx for idx, seq_id in enumerate(sequence_ids)}
    
    # Build sparse matrix by reading TSV in chunks to avoid memory issues
    row_indices = []
    col_indices = []
    distances = []
    
    # Check if file is empty first
    import os
    if os.path.getsize(results_tsv_path) == 0:
        print("Warning: Results file is empty, creating empty sparse matrix")
        sparse_matrix = coo_matrix((n_sequences, n_sequences))
    else:
        # Read header to check if file has content
        with open(results_tsv_path, 'r') as f:
            header = f.readline().strip()
            if not f.readline():  # Check if there's any data after header
                print("Warning: Results file has no data rows, creating empty sparse matrix")
                sparse_matrix = coo_matrix((n_sequences, n_sequences))
            else:
                # Process file in chunks
                print(f"Reading results in chunks of {chunk_size} rows")
                chunk_count = 0
                
                for chunk_df in pd.read_csv(results_tsv_path, sep="\t", chunksize=chunk_size):
                    chunk_count += 1
                    if chunk_count % 10 == 0:
                        print(f"Processing chunk {chunk_count}, total edges so far: {len(row_indices)}")
                    
                    # Convert ANI to distance
                    chunk_df['distance'] = (100 - chunk_df[score_type] * 100)
                    
                    # Map IDs to indices
                    chunk_df['i'] = chunk_df['query'].map(id_to_idx)
                    chunk_df['j'] = chunk_df['reference'].map(id_to_idx)
                    
                    # Remove any rows with unmapped IDs and self-comparisons
                    valid_mask = chunk_df['i'].notna() & chunk_df['j'].notna() & (chunk_df['i'] != chunk_df['j'])
                    chunk_df = chunk_df[valid_mask]
                    
                    if len(chunk_df) > 0:
                        # Append to lists
                        row_indices.extend(chunk_df['i'].astype(int).tolist())
                        col_indices.extend(chunk_df['j'].astype(int).tolist())
                        distances.extend(chunk_df['distance'].tolist())
                        
                        # Also add symmetric entries
                        row_indices.extend(chunk_df['j'].astype(int).tolist())
                        col_indices.extend(chunk_df['i'].astype(int).tolist())
                        distances.extend(chunk_df['distance'].tolist())
                
                print(f"Total edges loaded: {len(row_indices) // 2}")
                
                # Create sparse matrix
                if len(row_indices) > 0:
                    sparse_matrix = coo_matrix(
                        (distances, (row_indices, col_indices)), 
                        shape=(n_sequences, n_sequences)
                    )
                else:
                    print("Warning: No valid comparisons found, creating empty sparse matrix")
                    sparse_matrix = coo_matrix((n_sequences, n_sequences))
    
    # UMAP can handle sparse matrices directly
    print(f"Running UMAP with n_neighbors={config.n_neighbors}, min_dist={config.min_dist}")
    config.metric = 'precomputed'
    reducer = UMAP(**config.to_dict())
    embedding = reducer.fit_transform(sparse_matrix)
    
    print("UMAP embedding complete")
    
    return UMAPResult(
        embedding=embedding,
        sequence_ids=sequence_ids,
        config=config
    )
