from dataclasses import dataclass
import time
from typing import Iterator, List, Literal, NamedTuple
from datetime import datetime

from Bio.SeqRecord import SeqRecord
import numpy
from pandas.core.frame import DataFrame

from document_paths import DocumentPaths


class LzaniSettings(NamedTuple):
    exec_path: str
    score_type: Literal["ani"] | Literal["gani"] | Literal["tani"]


class ParasailSettings(NamedTuple):
    process_count: int


class RunSettings(NamedTuple):
    doc_paths: DocumentPaths
    analysis_method: Literal["parasail"] | Literal["lzani"]
    fasta_path: str
    output_path: str
    cluster_method: str  # Literal[reorder_methods]
    lzani: LzaniSettings
    parasail: ParasailSettings


class WorkflowResult(NamedTuple):
    seq_dict: dict[str, str]
    ordered_ids: List[str]
    reordered_ids: List[str]
    min_score: float
    max_sequence_length: int
    distance_matrix: numpy.ndarray
    similarity_matrix: DataFrame
    reordered_ids: List[str]
    is_aa: bool
    warnings: List[str]
    error: str | None


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
