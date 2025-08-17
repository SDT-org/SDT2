from dataclasses import dataclass
import time
from typing import List, Literal, NamedTuple
from datetime import datetime

import numpy
from pandas.core.frame import DataFrame

from document_paths import DocumentPaths


class LzaniSettings(NamedTuple):
    exec_path: str
    score_type: Literal["ani"]  | Literal["tani"]
    aw: int | None = None
    am: int | None = None
    mal: int | None = None
    msl: int | None = None
    mrd: float | None = None
    mqd: float | None = None
    reg: int | None = None
    ar: float | None = None


class ParasailSettings(NamedTuple):
    process_count: int
    scoring_matrix: str | None = None # None for custom matrix selection
    open_penalty: int | None = None
    extend_penalty: int | None = None


class VclustSettings(NamedTuple):
    # Prefiltering settings
    kmer_min_similarity: float = 0.30  # --min-ident in prefilter
    kmer_min_kmers: int = 20  # --min-kmers in prefilter
    kmer_fraction: float = 0.5  # --kmers-fraction in prefilter
    kmer_size: int = 25  # -k in prefilter
    batch_size: int = 0  # --batch-size in prefilter, 0 means no batching
    max_seqs_per_query: int = 0  # --max-seqs in prefilter, 0 means no limit
    
    # Clustering settings
    cluster_algorithm: str = "cd-hit"  # --algorithm in cluster
    cluster_metric: str = "tani"  # --metric in cluster
    cluster_threshold: float = 0.70  # --tani/--gani/--ani in cluster (also referred to as cdhit_threshold in frontend)
    min_coverage: float = 0.0  # --qcov/--rcov in cluster
    min_length_ratio: float = 0.0  # --len_ratio in cluster
    max_alignments: int = 0  # --num_alns in cluster, 0 means no limit
    use_representatives: bool = True  # --out-repr in cluster
    
    # LZ-ANI advanced settings
    lzani_mal: int = 11  # --mal in align
    lzani_msl: int = 7  # --msl in align
    lzani_mrd: int = 40  # --mrd in align
    lzani_mqd: int = 40  # --mqd in align
    lzani_reg: int = 35  # --reg in align
    lzani_aw: int = 15  # --aw in align
    lzani_am: int = 7  # --am in align
    lzani_ar: int = 3  # --ar in align
    
    # Output settings
    verbosity_level: int = 1  # 0=errors only, 1=info, 2=debug
    threads: int = 0  # 0 means use all available


class RunSettings(NamedTuple):
    doc_paths: DocumentPaths
    analysis_method: Literal["parasail"] | Literal["lzani"] | Literal["vclust"]
    fasta_path: str
    output_path: str
    cluster_method: str  # Literal[reorder_methods]
    lzani: LzaniSettings
    parasail: ParasailSettings
    vclust: VclustSettings | None = None
    export_alignments: bool = False
    alignment_export_path: str = ""


class WorkflowResult(NamedTuple):
    seq_dict: dict[str, str]
    ordered_ids: List[str]
    reordered_ids: List[str]
    min_score: float
    max_sequence_length: int
    distance_matrix: numpy.ndarray
    similarity_matrix: DataFrame
    reordered_ids: List[str]
    warnings: List[str]
    error: str | None
    is_aa: bool | None


@dataclass
class WorkflowRun:
    result: WorkflowResult
    settings: RunSettings
    progress: int | None
    stage: str = "Initializing"
    analyze_start_time: datetime | None = None
    analyze_start_counter: float | None = None

    def set_stage(self, stage: str):
        self.stage = stage
        if self.analyze_start_time is None and stage == "Analyzing":
            self.analyze_start_time = datetime.now()
            self.analyze_start_counter = time.perf_counter()
        print(f"Stage: {stage}")

    def set_progress(self, progress: int | None):
        old_progress = self.progress
        self.progress = progress
        if progress and old_progress is not None and progress > old_progress:   
            print(f"Progress: {progress}%")

    def valid(self) -> bool:
        return not self.result.error
