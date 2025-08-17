from typing import Protocol, Callable
from multiprocessing.sharedctypes import Synchronized

from workflow.models import RunSettings, WorkflowResult
from . import parasail
from . import lzani
from .. import vclust_pipeline


class AnalysisJob(Protocol):
    def run(
        self,
        result: WorkflowResult,
        settings: RunSettings,
        set_progress: Callable[[int], None],
        canceled: Synchronized,
    ) -> WorkflowResult:
        return result


jobs: dict[str, AnalysisJob] = {
    "parasail": parasail,
    "lzani": lzani,
    "vclust": vclust_pipeline,
}
