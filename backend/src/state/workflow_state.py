from typing import Dict
from workflow.models import WorkflowResult, WorkflowRun

# Global workflow state
workflow_runs: Dict[str, WorkflowRun] = {}
parsed_workflow_results: Dict[str, WorkflowResult] = {}
