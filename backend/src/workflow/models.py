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
    kmer_min_similarity: float = 0.30
    kmer_min_kmers: int = 2
    kmer_fraction: float = 0.5
    cdhit_threshold: float = 0.70


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
