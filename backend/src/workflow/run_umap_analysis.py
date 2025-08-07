#!/usr/bin/env python3
import sys
import os
import glob

# Add parent directory to path so imports work
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from workflow.umap_analysis import batch_analysis, compare_clustering_methods


def main(input_dir: str, output_dir: str):
    # Required files:
    # - distance_matrix.csv
    # - cluster_*.csv (cluster_ward.csv, cluster_complete.csv, etc.)
    # - metadata.csv (optional)
    
    distance_matrix_path = os.path.join(input_dir, 'distance_matrix.csv')
    if not os.path.exists(distance_matrix_path):
        print(f"Error: distance_matrix.csv not found in {input_dir}")
        sys.exit(1)
    
    # Find all cluster files
    cluster_files = glob.glob(os.path.join(input_dir, 'cluster_*.csv'))
    if not cluster_files:
        print(f"Warning: No cluster_*.csv files found in {input_dir}")
    
    cluster_csv_paths = {}
    for cluster_file in cluster_files:
        method_name = os.path.basename(cluster_file).replace('cluster_', '').replace('.csv', '')
        cluster_csv_paths[method_name] = cluster_file
    
    metadata_csv_path = os.path.join(input_dir, 'metadata.csv')
    if not os.path.exists(metadata_csv_path):
        metadata_csv_path = None
        print("Note: No metadata.csv found, proceeding without metadata")
    
    # Run batch analysis
    print(f"Running UMAP analysis...")
    print(f"  Distance matrix: {distance_matrix_path}")
    print(f"  Cluster methods: {list(cluster_csv_paths.keys())}")
    print(f"  Metadata: {'Yes' if metadata_csv_path else 'No'}")
    print(f"  Output directory: {output_dir}")
    
    results = batch_analysis(
        distance_matrix_path=distance_matrix_path,
        output_dir=output_dir,
        cluster_csv_paths=cluster_csv_paths if cluster_csv_paths else None,
        metadata_csv_path=metadata_csv_path,
        n_neighbors_list=[5, 15, 30],
        min_dist_list=[0.1, 0.25, 0.5]
    )
    
    print(f"\nCompleted! Generated {len(results)} UMAP embeddings")
    
    # Also create clustering comparison
    if cluster_csv_paths:
        comparison_output = os.path.join(output_dir, 'clustering_comparison.json')
        compare_clustering_methods(
            distance_matrix_path=distance_matrix_path,
            cluster_results_dir=input_dir,
            output_path=comparison_output,
            n_neighbors=15,
            min_dist=0.1
        )
        print(f"Clustering comparison saved to {comparison_output}")


if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python run_umap_analysis.py <input_dir> <output_dir>")
        print("\nRequired files in input_dir:")
        print("  - distance_matrix.csv")
        print("  - cluster_*.csv files (e.g., cluster_ward.csv, cluster_complete.csv)")
        print("  - metadata.csv (optional, with columns: sequence_id, genus, family, etc.)")
        sys.exit(1)
    
    input_dir = sys.argv[1]
    output_dir = sys.argv[2]
    
    if not os.path.exists(input_dir):
        print(f"Error: Input directory '{input_dir}' does not exist")
        sys.exit(1)
    
    main(input_dir, output_dir)
