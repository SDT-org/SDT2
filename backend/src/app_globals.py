from typing import Dict
from workflow.models import WorkflowResult, WorkflowRun

app_window = None
canceled = None
workflow_runs: Dict[str, WorkflowRun] = {}
parsed_workflow_results: Dict[str, WorkflowResult] = {}


def initialize_app_globals(window):
    global app_window, canceled, workflow_runs, parsed_workflow_results

    app_window = window
    canceled = None
    workflow_runs = {}
    parsed_workflow_results = {}


def get_app_window():
    return assert_app_window()


def assert_app_window():
    if app_window is None:
        raise RuntimeError("Application window not initialized")
    return app_window


def set_canceled(cancel_event):
    global canceled
    canceled = cancel_event


def get_canceled():
    return canceled


def get_workflow_runs():
    return workflow_runs


def get_parsed_workflow_results():
    return parsed_workflow_results
