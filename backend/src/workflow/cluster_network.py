import numpy as np
import networkx as nx
from typing import Dict, List, Tuple, Optional
import json


def create_similarity_network(
    distance_matrix: np.ndarray,
    sequence_ids: List[str],
    similarity_threshold: float = 95.0,
    min_similarity_filter: float = 50.0  # Filter out noise below this threshold
) -> nx.Graph:
    G = nx.Graph()
    n = len(sequence_ids)
    
    for i, seq_id in enumerate(sequence_ids):
        G.add_node(seq_id, index=i)
    
    distance_threshold = 100 - similarity_threshold
    max_distance_filter = 100 - min_similarity_filter
    edge_count = 0
    filtered_count = 0
    
    for i in range(n):
        for j in range(i + 1, n):
            distance = distance_matrix[i, j]
            
            # First filter: skip if below minimum meaningful similarity
            if distance > max_distance_filter:
                filtered_count += 1
                continue
                
            # Second filter: only create edge if above user's threshold
            if distance <= distance_threshold:
                similarity = 100 - distance
                G.add_edge(
                    sequence_ids[i], 
                    sequence_ids[j], 
                    weight=similarity,
                    distance=distance
                )
                edge_count += 1
    
    total_possible = n * (n - 1) // 2
    print(f"Created network with {n} nodes and {edge_count} edges")
    print(f"Filtered out {filtered_count} edges below {min_similarity_filter}% similarity")
    print(f"User threshold: {similarity_threshold}% (kept {edge_count}/{total_possible - filtered_count} meaningful edges)")
    return G


def run_connected_components_clustering(G: nx.Graph) -> Dict[str, int]:
    clusters = {}
    for cluster_id, component in enumerate(nx.connected_components(G), 1):
        for node in component:
            clusters[node] = cluster_id
    
    for node in G.nodes():
        if node not in clusters:
            clusters[node] = 0
    
    return clusters


def run_louvain_clustering(G: nx.Graph, resolution: float = 1.0) -> Dict[str, int]:
    try:
        import community as community_louvain
        
        # Check if graph has nodes
        if len(G) == 0:
            print("Warning: Empty graph, no nodes to cluster")
            return {}
        
        # Check if graph has edges
        if len(G.edges()) == 0:
            print("Warning: No edges in graph, all nodes will be singleton clusters")
            return {node: i for i, node in enumerate(G.nodes(), 1)}
        
        partition = community_louvain.best_partition(
            G, 
            weight='weight',
            resolution=resolution
        )
        
        # Debug: check what type partition is
        print(f"Partition type: {type(partition)}, length: {len(partition) if hasattr(partition, '__len__') else 'N/A'}")
        
        if not isinstance(partition, dict):
            print(f"Error: Expected dict from best_partition, got {type(partition)}")
            return run_connected_components_clustering(G)
        
        clusters = {}
        cluster_mapping = {}
        next_cluster_id = 1
        
        for node, cluster in partition.items():
            if cluster not in cluster_mapping:
                cluster_mapping[cluster] = next_cluster_id
                next_cluster_id += 1
            clusters[node] = cluster_mapping[cluster]
        
        print(f"Louvain clustering found {len(cluster_mapping)} clusters")
        return clusters
        
    except ImportError:
        print("Warning: python-louvain not installed, falling back to connected components")
        return run_connected_components_clustering(G)
    except Exception as e:
        print(f"Error in Louvain clustering: {type(e).__name__}: {e}")
        print("Falling back to connected components")
        return run_connected_components_clustering(G)


def compute_network_layout(
    G: nx.Graph,
    layout_method: str = "spring",
    iterations: int = 50
) -> Dict[str, Tuple[float, float]]:
    if layout_method == "spring":
        pos = nx.spring_layout(
            G, 
            k=1/np.sqrt(len(G.nodes())), 
            iterations=iterations,
            weight='weight'
        )
    elif layout_method == "kamada_kawai":
        if len(G) > 1000:
            print("Warning: Kamada-Kawai is slow for large networks, using spring layout")
            pos = nx.spring_layout(G, k=1/np.sqrt(len(G.nodes())), iterations=iterations)
        else:
            pos = nx.kamada_kawai_layout(G, weight='distance')
    elif layout_method == "forceatlas2":
        try:
            from fa2 import ForceAtlas2
            forceatlas2 = ForceAtlas2(
                outboundAttractionDistribution=True,
                linLogMode=False,
                adjustSizes=False,
                edgeWeightInfluence=1.0,
                jitterTolerance=1.0,
                barnesHutOptimize=True,
                barnesHutTheta=1.2,
                multiThreaded=0,
                scalingRatio=2.0,
                strongGravityMode=False,
                gravity=1.0,
                verbose=False
            )
            pos = forceatlas2.forceatlas2_networkx_layout(G, pos=None, iterations=iterations)
        except ImportError:
            print("Warning: ForceAtlas2 not installed, using spring layout")
            pos = nx.spring_layout(G, k=1/np.sqrt(len(G.nodes())), iterations=iterations)
    else:
        pos = nx.spring_layout(G, k=1/np.sqrt(len(G.nodes())), iterations=iterations)
    
    return dict(pos)  # Convert to regular dict to fix type issue


def get_network_data(
    distance_matrix: np.ndarray,
    sequence_ids: List[str],
    similarity_threshold: float = 95.0,
    clustering_method: str = "louvain",
    layout_method: str = "spring",
    resolution: float = 1.0,
    min_similarity_filter: float = 50.0
) -> Dict:
    G = create_similarity_network(
        distance_matrix, 
        sequence_ids, 
        similarity_threshold,
        min_similarity_filter=min_similarity_filter
    )
    
    if clustering_method == "louvain":
        clusters = run_louvain_clustering(G, resolution)
    else:
        clusters = run_connected_components_clustering(G)
    
    pos = compute_network_layout(G, layout_method)
    
    nodes = []
    for node in G.nodes():
        x, y = pos[node]
        nodes.append({
            "id": node,
            "x": float(x),
            "y": float(y),
            "cluster": clusters.get(node, 0),
            "degree": G.degree(node)  # type: ignore
        })
    
    edges = []
    for source, target, data in G.edges(data=True):
        edges.append({
            "source": source,
            "target": target,
            "weight": float(data.get('weight', 1.0))
        })
    
    cluster_counts = {}
    for cluster_id in clusters.values():
        cluster_counts[cluster_id] = cluster_counts.get(cluster_id, 0) + 1
    
    noise_count = cluster_counts.get(0, 0)
    actual_clusters = {k: v for k, v in cluster_counts.items() if k != 0}
    
    stats = {
        "total_nodes": len(G.nodes()),
        "total_edges": len(G.edges()),
        "total_clusters": len(actual_clusters),
        "noise_points": noise_count,
        "largest_cluster_size": max(actual_clusters.values()) if actual_clusters else 0,
        "smallest_cluster_size": min(actual_clusters.values()) if actual_clusters else 0,
        "average_degree": sum(dict(G.degree()).values()) / len(G) if len(G) > 0 else 0,  # type: ignore
        "density": nx.density(G),
        "cluster_sizes": cluster_counts
    }
    
    return {
        "nodes": nodes,
        "edges": edges,
        "stats": stats,
        "bounds": {
            "x": [min(n["x"] for n in nodes), max(n["x"] for n in nodes)],
            "y": [min(n["y"] for n in nodes), max(n["y"] for n in nodes)]
        }
    }
