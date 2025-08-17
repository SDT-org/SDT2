import os
import subprocess
import tempfile
from pathlib import Path
from typing import List, Tuple, Optional, NamedTuple
import logging

logger = logging.getLogger(__name__)


class KmerDbResult(NamedTuple):
    filtered_pairs: List[Tuple[str, str]]
    similarity_values: List[float]
    output_file: Optional[Path]


def run_kmerdb_prefilter(
    fasta_path: Path,
    output_path: Optional[Path] = None,
    min_similarity: float = 0.7,
    min_kmers: int = 20,
    kmer_fraction: float = 1.0,
    threads: Optional[int] = None
) -> KmerDbResult:
    if not fasta_path.exists():
        raise FileNotFoundError(f"Input file not found: {fasta_path}")
    
    kmerdb_exe = Path(__file__).parent.parent.parent.parent / "bin" / "kmer-db.exe"
    if not kmerdb_exe.exists():
        raise FileNotFoundError(f"kmer-db executable not found at {kmerdb_exe}")
    
    if threads is None:
        threads = os.cpu_count() or 1
    
    with tempfile.TemporaryDirectory() as temp_dir:
        temp_path = Path(temp_dir)
        
        # Step 1: Build k-mer database
        db_path = temp_path / "sequences.kdb"
        build_cmd = [
            str(kmerdb_exe), "build",
            "-multisample-fasta", str(fasta_path),
            str(db_path),
            "-t", str(threads),
            "-k", "25"
        ]
        if kmer_fraction < 1.0:
            build_cmd.extend(["-f", str(kmer_fraction)])
        
        logger.info(f"Running kmer-db build: {' '.join(build_cmd)}")
        result = subprocess.run(build_cmd, capture_output=True, text=True)
        if result.returncode != 0:
            raise RuntimeError(f"kmer-db build failed: {result.stderr}")

        # Step 2: Run all2all-sp to get sparse results
        all2all_path = temp_path / "all2all.txt"
        all2all_cmd = [
            str(kmerdb_exe), "all2all-sp",
            str(db_path),
            str(all2all_path),
            "-t", str(threads),
            "-min", f"num-kmers:{min_kmers}",
            "-min", f"ani-shorter:{min_similarity}",
            "-sparse"
        ]
        
        logger.info(f"Running kmer-db all2all-sp: {' '.join(all2all_cmd)}")
        result = subprocess.run(all2all_cmd, capture_output=True, text=True)
        if result.returncode != 0:
            raise RuntimeError(f"kmer-db all2all-sp failed: {result.stderr}")

        # Step 3: Run distance to get the final filter file
        distance_path = temp_path / "distance.txt"
        if output_path:
            distance_path = output_path
        
        distance_cmd = [
            str(kmerdb_exe), "distance", "ani-shorter",
            "-sparse",
            "-min", str(min_similarity),
            "-t", str(threads),
            str(all2all_path),
            str(distance_path)
        ]

        logger.info(f"Running kmer-db distance: {' '.join(distance_cmd)}")
        result = subprocess.run(distance_cmd, capture_output=True, text=True)
        if result.returncode != 0:
            raise RuntimeError(f"kmer-db distance failed: {result.stderr}")

        # Count the number of pairs in the final output file
        pair_count = 0
        with open(distance_path, 'r') as f:
            for line in f:
                if line.strip():
                    pair_count += 1
        
        logger.info(f"Prefiltering complete: {pair_count} pairs passed threshold")

        return KmerDbResult(
            filtered_pairs=[],  # This is no longer used, but kept for compatibility
            similarity_values=[],
            output_file=distance_path
        )


def estimate_memory_usage(
    num_sequences: int,
    avg_sequence_length: int,
    kmer_fraction: float = 1.0
) -> float:
    # Based on Vclust paper benchmarks
    base_memory_per_seq = 0.001  # GB
    kmer_memory_per_seq = (avg_sequence_length / 1000) * 0.0001 * kmer_fraction
    total_memory_gb = num_sequences * (base_memory_per_seq + kmer_memory_per_seq)
    
    # Sparse matrix overhead
    overhead_factor = 1.2 if num_sequences > 10000 else 1.5
    return total_memory_gb * overhead_factor
