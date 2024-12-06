import React from "react";
import type { AppState } from "../appState";

export const useStartProcessData = (appState: AppState) =>
  React.useCallback(() => {
    window.pywebview.api
      .run_process_data({
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
          window.pywebview.api.cancel_run();
        } else {
          throw e;
        }
      });
  }, [appState]);
