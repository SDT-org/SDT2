import json
import os
import numpy as np
import pandas as pd
from pandas import DataFrame
from numpy import eye, where, nan, nanmin, nanmax

from document_paths import build_document_paths
from file_utils import file_exists, read_json_file
from transformations import to_triangle, read_csv_matrix, read_stats_csv
from workflow import cluster
from app_state import get_document, update_document


class Data:
    def __init__(self):
        # Cache for UMAP distance matrices to avoid reloading on parameter changes
        self._umap_cache = {}
    
    def get_data(self, doc_id: str):
        doc = get_document(doc_id)
        doc_paths = build_document_paths(doc.tempdir_path)

        # Check if this is a large dataset that skipped matrix building
        is_large_dataset = False
        ids = []
        data = np.zeros((1, 1))  # Default to dummy matrix
        
        if file_exists(doc_paths.lzani_results_ids):
            # Check sequence count from lzani IDs file
            ids_df = pd.read_csv(doc_paths.lzani_results_ids, sep="\t")
            sequence_count = len(ids_df)
            if sequence_count > 2500:
                is_large_dataset = True
                ids = ids_df["id"].tolist()
                update_document(doc_id, sequences_count=sequence_count)
                # For large datasets, keep dummy matrix for UI compatibility
                # Real data will be loaded on-demand for UMAP
        
        if not is_large_dataset:
            try:
                df = read_csv_matrix(doc_paths.matrix)
                if not isinstance(df, DataFrame) or df.empty:
                    raise Exception(
                        f"Matrix file is empty or not a valid DataFrame: {doc_paths.matrix}"
                    )
                ids = df.index.tolist()
                update_document(doc_id, sequences_count=len(ids))
                data = df.to_numpy()
            except Exception as e:
                # If matrix reading fails, return minimal response
                print(f"Warning: Could not read matrix file: {e}")
                return json.dumps({
                    "metadata": {"minVal": 0, "maxVal": 100},
                    "data": [[], [[]]],  # Minimal valid data
                    "ids": [],
                    "identity_scores": [],
                    "full_stats": [],
                })

        if os.path.exists(doc_paths.stats):
            stats_df = read_stats_csv(doc_paths.stats)
        else:
            stats_df = DataFrame([])

        identity_scores = []
        all_scores = []
        id_map = {id: idx for idx, id in enumerate(ids)}

        is_lzani = False
        if file_exists(doc_paths.run_settings):
            run_settings = read_json_file(doc_paths.run_settings)
            if run_settings:
                is_lzani = run_settings.get("analysis_method") == "lzani"

        # For large datasets, skip expensive matrix operations
        if is_large_dataset:
            # Return minimal data for UI compatibility
            min_val = 0
            max_val = 100
            parsedData = [[]]
            unaligned_count = 0
        else:
            # Normal processing for smaller datasets
            for i in range(len(ids)):
                for j in range(i):
                    if i < data.shape[0] and j < data.shape[1]:
                        identity_score = 100 - data[i, j]
                        all_scores.append(identity_score)

                        if not (is_lzani and identity_score < 1):
                            identity_scores.append(
                                [id_map[ids[i]], id_map[ids[j]], identity_score]
                            )

            unaligned_count = len([s for s in all_scores if s < 3]) if is_lzani else 0

            heat_data = DataFrame(data, index=ids if ids else ["dummy"])
            heat_data = to_triangle(heat_data)

            if not isinstance(heat_data, DataFrame):
                raise Exception("Heat data is not a valid DataFrame.")

            heat_data_np = heat_data.to_numpy()
            diag_mask = eye(heat_data_np.shape[0], dtype=bool)
            heat_data_no_diag = where(diag_mask, nan, heat_data_np)

            if is_lzani:
                valid_values = heat_data_no_diag[
                    ~np.isnan(heat_data_no_diag) & (heat_data_no_diag >= 10)
                ]
                if valid_values.size > 0:
                    min_val = int(nanmin(valid_values))
                else:
                    min_val = int(nanmin(heat_data_no_diag)) if heat_data_no_diag.size > 0 else 0
            else:
                min_val = int(nanmin(heat_data_no_diag)) if heat_data_no_diag.size > 0 else 0

            max_val = int(nanmax(heat_data_no_diag)) if heat_data_no_diag.size > 0 else 100
            parsedData = heat_data.values.tolist()

        metadata = dict(minVal=min_val, maxVal=max_val)
        if file_exists(doc_paths.run_settings):
            metadata["run"] = read_json_file(doc_paths.run_settings)

        if is_lzani:
            metadata["unaligned_count"] = unaligned_count

        def calculate_stats(values):
            if len(values) == 0:
                return None
            return {
                "mean": float(np.mean(values)),
                "median": float(np.median(values)),
                "std": float(np.std(values)),
                "min": float(np.min(values)),
                "max": float(np.max(values)),
                "q1": float(np.percentile(values, 25)),
                "q3": float(np.percentile(values, 75)),
                "count": len(values),
            }

        if identity_scores:
            scores_array = np.array([score[2] for score in identity_scores])
            metadata["distribution_stats"] = calculate_stats(scores_array)

        if not stats_df.empty:
            stats_array = np.array(stats_df.values)
            metadata["gc_stats"] = calculate_stats(stats_array[:, 1])
            metadata["length_stats"] = calculate_stats(stats_array[:, 2])

        data_to_dump = dict(
            metadata=metadata,
            data=[ids] + parsedData,
            ids=ids,
            identity_scores=identity_scores,
            full_stats=stats_df.values.tolist(),
        )
        return json.dumps(data_to_dump)

    def get_clustermap_data(self, doc_id: str, threshold: float, method: str):
        doc = get_document(doc_id)

        doc_paths = build_document_paths(doc.tempdir_path)
        matrix_df = read_csv_matrix(doc_paths.matrix)
        matrix_np = matrix_df.to_numpy()
        matrix_np = np.round(matrix_np, 2)
        sorted_ids = matrix_df.index.tolist()

        seqid_clusters_df = cluster.get_clusters_dataframe(
            matrix_np, method, threshold, sorted_ids
        )
        seqid_clusters_df = seqid_clusters_df.rename(
            columns={
                str(seqid_clusters_df.columns[0]): "id",
                str(seqid_clusters_df.columns[1]): "cluster",
            }
        )

        seqid_clusters_df["id"] = seqid_clusters_df["id"].astype(str)
        seqid_clusters_df["original_cluster"] = seqid_clusters_df["cluster"]

        new_order = cluster.get_linkage_method_order(
            matrix_np, method, sorted_ids, threshold
        )

        reorder_indices = [sorted_ids.index(id) for id in new_order]

        reordered_matrix = matrix_np[np.ix_(reorder_indices, reorder_indices)]

        reordered_tick_text = [str(sorted_ids[i]) for i in reorder_indices]

        id_to_cluster = {
            row["id"]: row["cluster"] for _, row in seqid_clusters_df.iterrows()
        }

        seen_clusters = {}
        next_cluster_num = 1

        total_clusters = len(seqid_clusters_df["cluster"].unique())

        largest_cluster = max(seqid_clusters_df["cluster"].value_counts())

        cluster_counts = seqid_clusters_df["cluster"].value_counts()
        singleton_clusters = cluster_counts[cluster_counts == 1]
        singletons = len(singleton_clusters)
        print(singletons, "singletons")

        for seq_id in reordered_tick_text:
            original_cluster = id_to_cluster.get(seq_id)
            if original_cluster is not None and original_cluster not in seen_clusters:
                seen_clusters[original_cluster] = next_cluster_num
                next_cluster_num += 1

        for idx, row in seqid_clusters_df.iterrows():
            original_cluster = row["cluster"]
            if original_cluster in seen_clusters:
                seqid_clusters_df.at[idx, "cluster"] = seen_clusters[original_cluster]

        cluster_data = seqid_clusters_df.to_dict(orient="records")

        cluster_stats = {
            "total_clusters": total_clusters,
            "largest_cluster": int(largest_cluster),
            "singleton_clusters": singletons,
        }

        return {
            "matrix": to_triangle(reordered_matrix, fill_value=None).tolist(),
            "tickText": reordered_tick_text,
            "clusterData": cluster_data,
            "cluster_stats": cluster_stats,
        }

    def get_umap_data(self, doc_id: str, params: dict):
        from workflow.umap import run_umap, UMAPConfig
        from workflow.cluster_hdbscan import run_hdbscan_clustering, get_cluster_stats
        
        doc = get_document(doc_id)
        doc_paths = build_document_paths(doc.tempdir_path)
        
        # Check if we have cached distance matrix for this document
        if doc_id in self._umap_cache:
            print(f"Using cached distance matrix for doc {doc_id}")
            distance_matrix, sequence_ids = self._umap_cache[doc_id]
        else:
            print(f"Loading distance matrix for doc {doc_id}")
            # Check if this is lzani data and load directly from TSV for better performance
            if file_exists(doc_paths.run_settings):
                run_settings = read_json_file(doc_paths.run_settings)
                if run_settings and run_settings.get("analysis_method") == "lzani" and file_exists(doc_paths.lzani_results):
                    # Load lzani TSV directly without building full matrix
                    distance_matrix, sequence_ids = self._load_lzani_tsv_for_umap(
                        doc_paths.lzani_results,
                        doc_paths.lzani_results_ids,
                        run_settings.get("lzani", {}).get("score_type", "ani")
                    )
                else:
                    # Load from pre-computed matrix (parasail or existing lzani matrix)
                    matrix_df = read_csv_matrix(doc_paths.matrix)
                    distance_matrix = matrix_df.to_numpy()
                    sequence_ids = matrix_df.index.to_list()
            else:
                # Default: Load from matrix
                matrix_df = read_csv_matrix(doc_paths.matrix)
                distance_matrix = matrix_df.to_numpy()
                sequence_ids = matrix_df.index.to_list()
            
            # Cache the distance matrix for future use
            self._umap_cache[doc_id] = (distance_matrix, sequence_ids)
            print(f"Cached distance matrix for doc {doc_id} ({len(sequence_ids)} sequences)")
        
        # Run UMAP with specified parameters
        config = UMAPConfig(
            n_neighbors=params.get("n_neighbors", 15),
            min_dist=params.get("min_dist", 0.1),
            metric='precomputed'
        )
        
        umap_result = run_umap(
            distance_matrix=distance_matrix,
            sequence_ids=sequence_ids,
            config=config
        )
        
        # Run HDBSCAN clustering
        # Support both old and new parameter names for backward compatibility
        min_cluster_size = params.get("min_cluster_size", params.get("threshold", 5))
        cluster_epsilon = params.get("cluster_epsilon", 0.0)
        
        # Extract epsilon from methods array if using old format
        methods = params.get("methods", [])
        if methods and len(methods) > 0 and methods[0].startswith("hdbscan-"):
            try:
                # Epsilon comes as similarity percentage, convert to distance
                similarity_epsilon = float(methods[0].split("-")[1])
                cluster_epsilon = 100 - similarity_epsilon  # Convert similarity to distance
            except:
                pass
        
        cluster_assignments = run_hdbscan_clustering(
            distance_matrix=distance_matrix,
            sequence_ids=sequence_ids,
            min_cluster_size=min_cluster_size,
            cluster_selection_epsilon=cluster_epsilon
        )
        
        # Get cluster statistics
        cluster_stats = get_cluster_stats(cluster_assignments)
        
        # Debug logging
        print(f"HDBSCAN clustering: min_cluster_size={min_cluster_size}, epsilon={cluster_epsilon}")
        print(f"Distance matrix shape: {distance_matrix.shape}")
        print(f"Distance matrix range: min={distance_matrix.min():.3f}, max={distance_matrix.max():.3f}")
        print(f"Total sequences: {cluster_stats['total_sequences']}")
        print(f"Total clusters: {cluster_stats['total_clusters']}")
        print(f"Noise points: {cluster_stats['noise_points']}")
        
        # Format data for frontend
        embedding_data = []
        for i, seq_id in enumerate(sequence_ids):
            embedding_data.append({
                "id": seq_id,
                "x": float(umap_result.embedding[i, 0]),
                "y": float(umap_result.embedding[i, 1]),
                "cluster": int(cluster_assignments.get(seq_id, 0)),
                "clusters": {"hdbscan": int(cluster_assignments.get(seq_id, 0))}  # Keep compatibility
            })
        
        # Calculate bounds with padding
        x_coords = umap_result.embedding[:, 0]
        y_coords = umap_result.embedding[:, 1]
        x_range = x_coords.max() - x_coords.min()
        y_range = y_coords.max() - y_coords.min()
        padding = 0.1  # 10% padding
        
        return {
            "data": {
                "embedding": embedding_data,
                "bounds": {
                    "x": [float(x_coords.min() - x_range * padding), float(x_coords.max() + x_range * padding)],
                    "y": [float(y_coords.min() - y_range * padding), float(y_coords.max() + y_range * padding)]
                },
                "clusterStats": cluster_stats,
                "min_cluster_size": int(min_cluster_size),
                "cluster_epsilon": float(cluster_epsilon)
            },
            "metadata": {
                "n_neighbors": int(config.n_neighbors),
                "min_dist": float(config.min_dist),
                "sequences_count": int(len(sequence_ids)),
                "cluster_method": "hdbscan",
                "min_cluster_size": int(min_cluster_size),
                "cluster_epsilon": float(cluster_epsilon)
            }
        }

    def upload_metadata(self, doc_id: str, csv_content: str):
        """
        Upload and process metadata CSV for UMAP visualization.
        Returns available columns and match statistics.
        """
        doc = get_document(doc_id)
        doc_paths = build_document_paths(doc.tempdir_path)
        
        # Parse CSV
        try:
            import io
            metadata_df = pd.read_csv(io.StringIO(csv_content))
        except Exception as e:
            raise Exception(f"Invalid CSV format: {str(e)}")
        
        if len(metadata_df.columns) < 2:
            raise Exception("CSV must have at least 2 columns (ID column + metadata)")
        
        # First column is assumed to be sequence IDs
        id_column = metadata_df.columns[0]
        metadata_ids = metadata_df[id_column].astype(str).tolist()
        
        # Load sequence IDs - for large datasets, get from lzani IDs file
        if file_exists(doc_paths.lzani_results_ids):
            # Large dataset - get IDs from lzani results
            ids_df = pd.read_csv(doc_paths.lzani_results_ids, sep="\t")
            sequence_ids = ids_df["id"].tolist()
        else:
            # Regular dataset - get from matrix
            matrix_df = read_csv_matrix(doc_paths.matrix)
            sequence_ids = matrix_df.index.tolist()
        
        # Perform ID matching
        matches, match_stats = self._match_sequence_ids(metadata_ids, sequence_ids)
        
        # Save metadata to temp directory
        metadata_path = os.path.join(doc.tempdir_path, "metadata.csv")
        metadata_df.to_csv(metadata_path, index=False)
        
        # Save matching information
        match_info_path = os.path.join(doc.tempdir_path, "metadata_matches.json")
        with open(match_info_path, 'w') as f:
            json.dump({
                "matches": matches,
                "stats": match_stats,
                "id_column": id_column,
                "columns": metadata_df.columns.tolist()[1:]  # Exclude ID column
            }, f)
        
        # Detect column types
        column_info = {}
        for col in metadata_df.columns[1:]:  # Skip ID column
            # Simple type detection
            try:
                pd.to_numeric(metadata_df[col])
                column_info[col] = "numeric"
            except:
                column_info[col] = "categorical"
        
        return {
            "columns": metadata_df.columns.tolist()[1:],
            "column_types": column_info,
            "match_stats": match_stats,
            "total_rows": len(metadata_df)
        }
    
    def _match_sequence_ids(self, metadata_ids, sequence_ids):
        """
        Match metadata IDs to sequence IDs with GenBank version handling.
        Returns: (matches dict, stats dict)
        """
        matches = {}
        exact_matches = 0
        version_matches = 0
        unmatched = 0
        
        # Create lookup sets for faster matching
        sequence_id_set = set(sequence_ids)
        sequence_base_ids = {}
        
        # Build base ID map (without version numbers)
        for seq_id in sequence_ids:
            if '.' in seq_id:
                base_id = seq_id.split('.')[0]
                sequence_base_ids[base_id] = seq_id
        
        # Match metadata IDs
        for meta_id in metadata_ids:
            if meta_id in sequence_id_set:
                # Exact match
                matches[meta_id] = meta_id
                exact_matches += 1
            elif '.' in meta_id:
                # Try matching without version
                base_id = meta_id.split('.')[0]
                if base_id in sequence_base_ids:
                    matches[meta_id] = sequence_base_ids[base_id]
                    version_matches += 1
                else:
                    unmatched += 1
            else:
                # Try to find versioned match
                if meta_id in sequence_base_ids:
                    matches[meta_id] = sequence_base_ids[meta_id]
                    version_matches += 1
                else:
                    unmatched += 1
        
        stats = {
            "total_metadata_ids": len(metadata_ids),
            "total_sequence_ids": len(sequence_ids),
            "exact_matches": exact_matches,
            "version_matches": version_matches,
            "unmatched": unmatched,
            "match_percentage": round((exact_matches + version_matches) / len(sequence_ids) * 100, 1) if len(sequence_ids) > 0 else 0
        }
        
        return matches, stats
    
    def get_metadata_for_umap(self, doc_id: str, column_name: str):

        doc = get_document(doc_id)
        
        # Load metadata and matching info
        metadata_path = os.path.join(doc.tempdir_path, "metadata.csv")
        match_info_path = os.path.join(doc.tempdir_path, "metadata_matches.json")
        
        if not os.path.exists(metadata_path) or not os.path.exists(match_info_path):
            raise Exception("No metadata uploaded")
        
        metadata_df = pd.read_csv(metadata_path)
        with open(match_info_path, 'r') as f:
            match_info = json.load(f)
        
        # Get the requested column
        id_column = match_info["id_column"]
        if column_name not in metadata_df.columns:
            raise Exception(f"Column '{column_name}' not found in metadata")
        
        # Create value map using the matching information
        value_map = {}
        for _, row in metadata_df.iterrows():
            meta_id = str(row[id_column])
            if meta_id in match_info["matches"]:
                seq_id = match_info["matches"][meta_id]
                value_map[seq_id] = row[column_name]
        
        return {
            "column_name": column_name,
            "value_map": value_map,
            "column_type": "numeric" if pd.api.types.is_numeric_dtype(metadata_df[column_name]) else "categorical"
        }
    
    def _load_lzani_tsv_for_umap(self, results_tsv_path: str, ids_tsv_path: str, score_column: str = "ani"):

        ids_df = pd.read_csv(ids_tsv_path, sep="\t")
        all_ids = ids_df["id"].tolist()
        n_sequences = len(all_ids)
        
        results_df = pd.read_csv(results_tsv_path, sep="\t")
        
        # Create pivot table with similarity scores
        if results_df.empty:
            matrix = pd.DataFrame(index=all_ids, columns=all_ids, data=100.0)
        else:
            matrix = results_df.pivot_table(
                index="query", columns="reference", values=score_column, aggfunc="first"
            )
            # Reindex to ensure all IDs are present
            matrix = matrix.reindex(index=all_ids, columns=all_ids)
        
        # Convert to numpy and scale to percentage
        matrix_np = matrix.to_numpy() * 100
        
        # Handle NaN values and convert similarity to distance
        matrix_np = np.where(np.isnan(matrix_np), 100, 100 - matrix_np)
        
        # Make symmetric (average forward and reverse scores)
        matrix_np = (matrix_np + matrix_np.T) / 2
        
        # Set diagonal to 0 (self-distance)
        np.fill_diagonal(matrix_np, 0)
        
        print(f"Loaded lzani TSV directly: {n_sequences} sequences")
        
        return matrix_np, all_ids
