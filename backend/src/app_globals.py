from typing import Dict
from workflow.models import WorkflowResult, WorkflowRun

app_window = None
js_api = None
canceled = None
workflow_runs: Dict[str, WorkflowRun] = {}
parsed_workflow_results: Dict[str, WorkflowResult] = {}


def initialize_app_globals(window, api):
    global app_window, canceled, workflow_runs, parsed_workflow_results, js_api

    app_window = window
    js_api = api
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
