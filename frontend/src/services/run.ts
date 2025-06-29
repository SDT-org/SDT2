import type { AppState } from "../appState";
import messages from "../messages";

export const startRun = (docId: string, appState: AppState) =>
  window.pywebview.api
    .start_workflow_run({
      doc_id: docId,
      cluster_method: appState.enableClustering
        ? appState.cluster_method
        : "None",
      compute_cores: appState.compute_cores,
      analysisMethod: appState.analysisMethod,
    })
    .catch((e) => {
      if (e.toString().includes("PARASAIL_TRACEBACK")) {
        alert(
          "An error occured while aligning. " +
            "Please ensure you have adequate swap/page size and system memory.\n\n" +
            "Error ID: PARASAIL_TRACEBACK",
        );
        cancelRun(docId, "preserve");
      } else if (
        e &&
        typeof e === "object" &&
        "message" in e &&
        typeof e.message === "string" &&
        Object.keys(messages).includes(e.message)
      ) {
        alert(`Error: ${messages[e.message as keyof typeof messages]}`);
        cancelRun(docId, "preserve");
      } else {
        throw e;
      }
    });

export const cancelRun = (
  docId: string,
  run_settings: Parameters<typeof window.pywebview.api.cancel_run>[1],
) => window.pywebview.api.cancel_run(docId, run_settings);
