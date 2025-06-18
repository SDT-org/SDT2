import type { AppState } from "../appState";

export const startRun = (docId: string, appState: AppState) =>
  window.pywebview.api
    .start_run({
      doc_id: docId,
      cluster_method: appState.enableClustering
        ? appState.cluster_method
        : "None",
      compute_cores: appState.compute_cores,
      export_alignments: appState.enableOutputAlignments ? "True" : "False",
      analysisMethod: appState.analysisMethod,
      lzani_score_type: appState.lzaniScoreType,
    })
    .catch((e) => {
      if (e.toString().includes("PARASAIL_TRACEBACK")) {
        alert(
          "An error occured while aligning. " +
            "Please ensure you have adequate swap/page size and system memory.\n\n" +
            "Error ID: PARASAIL_TRACEBACK",
        );
        cancelRun(docId, "preserve");
      } else {
        throw e;
      }
    });

export const cancelRun = (
  docId: string,
  run_settings: Parameters<typeof window.pywebview.api.cancel_run>[1],
) => window.pywebview.api.cancel_run(docId, run_settings);
