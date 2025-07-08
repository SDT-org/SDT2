import type { DocState } from "../appState";
import messages from "../messages";

export const startRun = (docState: DocState) =>
  window.pywebview.api
    .start_workflow_run({
      doc_id: docState.id,
      cluster_method: docState.enableClustering
        ? docState.cluster_method
        : "None",
      compute_cores: docState.compute_cores || 1,
      analysisMethod: docState.analysisMethod,
      ...(docState.analysisMethod === "parasail" && docState.overrideParasail
        ? docState.parasail_settings
        : {}),
      ...(docState.analysisMethod === "lzani" && docState.overrideLzani
        ? docState.lzani_settings
        : {}),
      lzani_score_type: docState.lzaniScoreType,
    })
    .catch((e) => {
      if (e.toString().includes("PARASAIL_TRACEBACK")) {
        alert(
          "An error occured while aligning. " +
            "Please ensure you have adequate swap/page size and system memory.\n\n" +
            "Error ID: PARASAIL_TRACEBACK",
        );
        cancelRun(docState.id, "preserve");
      } else if (
        e &&
        typeof e === "object" &&
        "message" in e &&
        typeof e.message === "string" &&
        Object.keys(messages).includes(e.message)
      ) {
        alert(`Error: ${messages[e.message as keyof typeof messages]}`);
        cancelRun(docState.id, "preserve");
      } else {
        throw e;
      }
    });

export const cancelRun = (
  docId: string,
  run_settings: Parameters<typeof window.pywebview.api.cancel_run>[1],
) => window.pywebview.api.cancel_run(docId, run_settings);
