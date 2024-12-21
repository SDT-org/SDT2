import type { AppState } from "../appState";

export const startRun = (appState: AppState) =>
  window.pywebview.api
    .start_run({
      cluster_method: appState.client.enableClustering
        ? appState.client.cluster_method
        : "None",
      compute_cores: appState.client.compute_cores,
      export_alignments: appState.client.enableOutputAlignments
        ? "True"
        : "False",
    })
    .catch((e) => {
      if (e.toString().includes("PARASAIL_TRACEBACK")) {
        alert(
          "An error occured while aligning. " +
            "Please ensure you have adequate swap/page size and system memory.\n\n" +
            "Error ID: PARASAIL_TRACEBACK",
        );
        cancelRun("preserve");
      } else {
        throw e;
      }
    });

export const cancelRun = (
  run_settings: Parameters<typeof window.pywebview.api.cancel_run>[0],
) => window.pywebview.api.cancel_run(run_settings);
