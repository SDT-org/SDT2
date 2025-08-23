import os
import sys
import json

current_file_path = os.path.abspath(os.path.join(os.path.dirname(__file__)))
sys.path.append(os.path.join(current_file_path, "../../"))
sys.path.append(os.path.join(current_file_path, "."))

from mime_setup import register_mimetypes

from app_state import (
    initialize_app_state,
    get_state,
    set_state,
    update_document,
)
from app_globals import (
    initialize_app_globals,
    get_workflow_runs,
    get_canceled,
)
import platform
import psutil
os.environ['QTWEBENGINE_REMOTE_DEBUGGING'] = '8228'
import webview
import json
from config import default_window_title, cpu_count
from api import system, files, documents, workflow, data, export, windows
from file_utils import get_html_path

register_mimetypes()


def do_cancel_run():
    canceled = get_canceled()
    if canceled:
        # It's ok if the manager has already released the resource
        try:
            canceled.value = True
        except:
            pass
    doc_id = get_state().active_run_document_id
    if doc_id:
        workflow_runs = get_workflow_runs()
        if doc_id in workflow_runs:
            del workflow_runs[doc_id]
        update_document(
            doc_id,
            view="runner",
            doc_id=doc_id,
            progress=0,
            pair_progress=0,
            pair_count=0,
        )
    set_state(active_run_document_id=None)
    print("Run canceled")


class Api:
    def __init__(self):
        self.system = system.System()
        self.files = files.Files()
        self.documents = documents.Documents()
        self.workflow = workflow.Workflow()
        self.data = data.Data()
        self.export = export.Export()
        self.windows = windows.Windows()


def push_backend_state(window: webview.Window):
    state_dict = get_state()._asdict()
    state_dict["documents"] = [d._asdict() for d in state_dict["documents"]]
    js_app_state = json.dumps({"state": state_dict})
    window.evaluate_js(
        f"document.dispatchEvent(new CustomEvent('sync-state', {{ detail: {js_app_state} }}))"
    )


def on_closed():
    do_cancel_run()
    os._exit(0)


if __name__ == "__main__":
    api = Api()
    app_window = webview.create_window(
        default_window_title,
        url=get_html_path(),
        js_api=api,
        # TODO: store last window size and position
        width=1200,
        height=900,
        min_size=(900, 870),
        confirm_close=False,
        # frameless=True,
        # easy_drag=False,
        # maximized=True
    )
    app_window.events.closed += on_closed

    initialize_app_globals(window=app_window, api=api)

    initialize_app_state(
        debug=os.getenv("DEBUG", "false").lower() == "true",
        platform=dict(
            platform=platform.platform(),
            cores=cpu_count,
            memory=psutil.virtual_memory().total,
        ),
        on_update=lambda _: push_backend_state(app_window),
    )
    webview.settings['OPEN_DEVTOOLS_IN_DEBUG'] = False
    webview.start(debug=get_state().debug, private_mode=False)
