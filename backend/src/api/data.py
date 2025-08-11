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
        self._umap_cache = {}
    
    def get_data(self, doc_id: str):
        doc = get_document(doc_id)
        doc_paths = build_document_paths(doc.tempdir_path)

        is_large_dataset = False
        ids = []
        data = np.zeros((1, 1))
        min_val = 0
        max_val = 100
        unaligned_count = 0
        parsedData = []
        
        if file_exists(doc_paths.lzani_results_ids):
            ids_df = pd.read_csv(doc_paths.lzani_results_ids, sep="\t")
            sequence_count = len(ids_df)
            if sequence_count > 2500:
                # Return error for large datasets - no heatmap/clustermap support
                update_document(doc_id, sequences_count=sequence_count)
                return json.dumps({
                    "error": "Dataset too large for heatmap visualization",
                    "sequences_count": sequence_count,
                    "message": "Please use UMAP or Network visualization for datasets with more than 2500 sequences"
                })
        
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
                print(f"Warning: Could not read matrix file: {e}")
                return json.dumps({
                    "metadata": {"minVal": 0, "maxVal": 100},
                    "data": [[], [[]]],
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

        if not is_large_dataset:
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
            dist_stats = calculate_stats(scores_array)
            if dist_stats:
                metadata["distribution_stats"] = dist_stats

        if not stats_df.empty:
            stats_array = np.array(stats_df.values)
            gc_stats = calculate_stats(stats_array[:, 1])
            length_stats = calculate_stats(stats_array[:, 2])
            if gc_stats:
                metadata["gc_stats"] = gc_stats
            if length_stats:
                metadata["length_stats"] = length_stats

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

        largest_cluster = int(seqid_clusters_df["cluster"].value_counts().max())

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
        
        if doc_id in self._umap_cache:
            print(f"Using cached distance matrix for doc {doc_id}")
            distance_matrix, sequence_ids = self._umap_cache[doc_id]
        else:
            print(f"Loading distance matrix for doc {doc_id}")
            if file_exists(doc_paths.run_settings):
                run_settings = read_json_file(doc_paths.run_settings)
                if run_settings and run_settings.get("analysis_method") == "lzani" and file_exists(doc_paths.lzani_results):
                    distance_matrix, sequence_ids = self._load_lzani_tsv_for_umap(
                        doc_paths.lzani_results,
                        doc_paths.lzani_results_ids,
                        run_settings.get("lzani", {}).get("score_type", "ani")
                    )
                else:
                    matrix_df = read_csv_matrix(doc_paths.matrix)
                    distance_matrix = matrix_df.to_numpy()
                    sequence_ids = matrix_df.index.to_list()
            else:
                matrix_df = read_csv_matrix(doc_paths.matrix)
                distance_matrix = matrix_df.to_numpy()
                sequence_ids = matrix_df.index.to_list()
            
            self._umap_cache[doc_id] = (distance_matrix, sequence_ids)
            print(f"Cached distance matrix for doc {doc_id} ({len(sequence_ids)} sequences)")
        
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
        
        min_cluster_size = params.get("min_cluster_size", params.get("threshold", 5))
        cluster_epsilon = params.get("cluster_epsilon", 0.0)
        
        methods = params.get("methods", [])
        if methods and len(methods) > 0 and methods[0].startswith("hdbscan-"):
            try:
                similarity_epsilon = float(methods[0].split("-")[1])
                cluster_epsilon = 100 - similarity_epsilon
            except:
                pass
        
        cluster_assignments = run_hdbscan_clustering(
            distance_matrix=distance_matrix,
            sequence_ids=sequence_ids,
            min_cluster_size=min_cluster_size,
            cluster_selection_epsilon=cluster_epsilon
        )
        
        cluster_stats = get_cluster_stats(cluster_assignments)
        
        print(f"HDBSCAN clustering: min_cluster_size={min_cluster_size}, epsilon={cluster_epsilon}")
        print(f"Distance matrix shape: {distance_matrix.shape}")
        print(f"Distance matrix range: min={distance_matrix.min():.3f}, max={distance_matrix.max():.3f}")
        print(f"Total sequences: {cluster_stats['total_sequences']}")
        print(f"Total clusters: {cluster_stats['total_clusters']}")
        print(f"Noise points: {cluster_stats['noise_points']}")
        
        embedding_data = []
        for i, seq_id in enumerate(sequence_ids):
            embedding_data.append({
                "id": seq_id,
                "x": float(umap_result.embedding[i, 0]),
                "y": float(umap_result.embedding[i, 1]),
                "cluster": int(cluster_assignments.get(seq_id, 0)),
                "clusters": {"hdbscan": int(cluster_assignments.get(seq_id, 0))}
            })
        
        x_coords = umap_result.embedding[:, 0]
        y_coords = umap_result.embedding[:, 1]
        x_range = x_coords.max() - x_coords.min()
        y_range = y_coords.max() - y_coords.min()
        padding = 0.1
        
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
        doc = get_document(doc_id)
        doc_paths = build_document_paths(doc.tempdir_path)
        
        try:
            import io
            metadata_df = pd.read_csv(io.StringIO(csv_content))
        except Exception as e:
            raise Exception(f"Invalid CSV format: {str(e)}")
        
        if len(metadata_df.columns) < 2:
            raise Exception("CSV must have at least 2 columns (ID column + metadata)")
        
        id_column = metadata_df.columns[0]
        metadata_ids = metadata_df[id_column].astype(str).tolist()
        
        if file_exists(doc_paths.lzani_results_ids):
            ids_df = pd.read_csv(doc_paths.lzani_results_ids, sep="\t")
            sequence_ids = ids_df["id"].tolist()
        else:
            matrix_df = read_csv_matrix(doc_paths.matrix)
            sequence_ids = matrix_df.index.tolist()
        
        matches, match_stats = self._match_sequence_ids(metadata_ids, sequence_ids)
        
        metadata_path = os.path.join(doc.tempdir_path, "metadata.csv")
        metadata_df.to_csv(metadata_path, index=False)
        
        match_info_path = os.path.join(doc.tempdir_path, "metadata_matches.json")
        with open(match_info_path, 'w') as f:
            json.dump({
                "matches": matches,
                "stats": match_stats,
                "id_column": id_column,
                "columns": metadata_df.columns.tolist()[1:]
            }, f)
        
        column_info = {}
        for col in metadata_df.columns[1:]:
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
        matches = {}
        exact_matches = 0
        version_matches = 0
        unmatched = 0
        
        sequence_id_set = set(sequence_ids)
        sequence_base_ids = {}
        
        for seq_id in sequence_ids:
            if '.' in seq_id:
                base_id = seq_id.split('.')[0]
                sequence_base_ids[base_id] = seq_id
        
        for meta_id in metadata_ids:
            if meta_id in sequence_id_set:
                matches[meta_id] = meta_id
                exact_matches += 1
            elif '.' in meta_id:
                base_id = meta_id.split('.')[0]
                if base_id in sequence_base_ids:
                    matches[meta_id] = sequence_base_ids[base_id]
                    version_matches += 1
                else:
                    unmatched += 1
            else:
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
    
    def get_network_data(self, doc_id: str, params: dict):
        from workflow.cluster_network import get_network_data
        
        doc = get_document(doc_id)
        doc_paths = build_document_paths(doc.tempdir_path)
        
        if file_exists(doc_paths.lzani_results_ids):
            run_settings = read_json_file(doc_paths.run_settings) if file_exists(doc_paths.run_settings) else {}
            score_type = run_settings.get("lzani", {}).get("score_type", "ani") if run_settings else "ani"
            distance_matrix, sequence_ids = self._load_lzani_tsv_for_umap(
                doc_paths.lzani_results,
                doc_paths.lzani_results_ids,
                score_type
            )
        else:
            matrix_df = read_csv_matrix(doc_paths.matrix)
            distance_matrix = matrix_df.to_numpy()
            sequence_ids = matrix_df.index.to_list()
        
        network_result = get_network_data(
            distance_matrix=distance_matrix,
            sequence_ids=sequence_ids,
            similarity_threshold=params.get("similarity_threshold", 95.0),
            clustering_method=params.get("clustering_method", "louvain"),
            layout_method=params.get("layout_method", "spring"),
            resolution=params.get("resolution", 1.0),
            min_similarity_filter=params.get("min_similarity_filter", 50.0)
        )
        
        return {
            "data": network_result,
            "metadata": {
                "similarity_threshold": params.get("similarity_threshold", 95.0),
                "clustering_method": params.get("clustering_method", "louvain"),
                "sequences_count": len(sequence_ids)
            }
        }
    
    def get_metadata_for_umap(self, doc_id: str, column_name: str):

        doc = get_document(doc_id)
        
        metadata_path = os.path.join(doc.tempdir_path, "metadata.csv")
        match_info_path = os.path.join(doc.tempdir_path, "metadata_matches.json")
        
        if not os.path.exists(metadata_path) or not os.path.exists(match_info_path):
            raise Exception("No metadata uploaded")
        
        metadata_df = pd.read_csv(metadata_path)
        with open(match_info_path, 'r') as f:
            match_info = json.load(f)
        
        id_column = match_info["id_column"]
        if column_name not in metadata_df.columns:
            raise Exception(f"Column '{column_name}' not found in metadata")
        
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
        
        if results_df.empty:
            matrix = pd.DataFrame(index=all_ids, columns=all_ids, data=100.0)
        else:
            matrix = results_df.pivot_table(
                index="query", columns="reference", values=score_column, aggfunc="first"
            )
            matrix = matrix.reindex(index=all_ids, columns=all_ids)
        
        matrix_np = matrix.to_numpy() * 100
        
        matrix_np = np.where(np.isnan(matrix_np), 100, 100 - matrix_np)
        
        # Take the maximum of forward and reverse scores
        matrix_np = np.maximum(matrix_np, matrix_np.T)
        
        np.fill_diagonal(matrix_np, 0)
        
        print(f"Loaded lzani TSV directly: {n_sequences} sequences")
        
        return matrix_np, all_ids
