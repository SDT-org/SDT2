"""
Cache for storing sparse matrices from vclust pipeline.
This avoids JSON serialization issues with scipy sparse matrices.
"""
from typing import Dict, Optional, Tuple
from scipy.sparse import coo_matrix

# Global cache for vclust sparse matrices
_vclust_cache: Dict[str, Tuple[coo_matrix, list]] = {}

def store_vclust_matrix(doc_id: str, distance_matrix: coo_matrix, ordered_ids: list):
    """Store a sparse matrix and ordered IDs for a document."""
    _vclust_cache[doc_id] = (distance_matrix, ordered_ids)

def get_vclust_matrix(doc_id: str) -> Optional[Tuple[coo_matrix, list]]:
    """Retrieve a sparse matrix and ordered IDs for a document."""
    return _vclust_cache.get(doc_id)

def clear_vclust_cache(doc_id: Optional[str] = None):
    """Clear the cache for a specific document or all documents."""
    if doc_id:
        _vclust_cache.pop(doc_id, None)
    else:
        _vclust_cache.clear()
