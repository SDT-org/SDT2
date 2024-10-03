import React, { ErrorInfo } from "react";
import messages from "./messages";

export const clusterMethods = ["Neighbor-Joining", "UPGMA", "None"] as const;

export type SaveableImageFormat = "png" | "jpeg" | "svg";

export type PerformanceProfile =
  | "recommended"
  | "best"
  | "balanced"
  | "low"
  | "custom";

export type AppState = {
  view: "runner" | "loader" | "viewer";
  filename: string[];
  filetype: string;
  basename: string;
  progress: number;
  debug: boolean;
  sequences_count: number;
  alignment_output_path: string;
  export_path: string;
  stage: string;
  pair_progress: number;
  pair_count: number;
  estimated_time?: number;
  validation_error_id?: keyof typeof messages;
  compute_stats?: {
    recommended_cores: number;
    required_memory: number;
    available_memory: number;
  };
  platform: {
    cores: number;
    memory: number;
    platform: string;
  };
  client: {
    dataView: "heatmap" | "plot";
    performanceProfile: PerformanceProfile;
    cluster_method: (typeof clusterMethods)[number];
    compute_cores: number;
    error?: Error | null;
    errorInfo?: ErrorInfo | PromiseRejectionEvent["reason"] | null;
    saveFormat: SaveableImageFormat;
    showExportModal: boolean;
  };
};

export const initialAppState: AppState = {
  view: "runner",
  filename: [""],
  filetype: "",
  basename: "",
  progress: 0,
  debug: false,
  sequences_count: 0,
  alignment_output_path: "",
  export_path: "",
  stage: "Preprocessing",
  pair_progress: 0,
  pair_count: 0,
  client: {
    dataView: "heatmap",
    cluster_method: "Neighbor-Joining",
    saveFormat: "svg",
    showExportModal: false,
    compute_cores: 1,
    performanceProfile: "recommended",
  },
  platform: {
    cores: 1,
    memory: 1,
    platform: "unknown",
  },
};

export type SetAppState = React.Dispatch<React.SetStateAction<AppState>>;

export const AppStateContext = React.createContext<{
  appState: AppState;
  setAppState: SetAppState;
}>({
  appState: initialAppState,
  setAppState: () => undefined,
});

export const useAppState = () => {
  const context = React.useContext(AppStateContext);
  if (context === undefined) {
    throw new Error("useAppState must be used within an AppStateProvider");
  }
  return context;
};

export const syncAppState = (setAppState: SetAppState) => {
  if (!window.pywebview) {
    console.warn("Frontend-only mode detected, app state will not be synced.");
    return Promise.resolve();
  }

  return window.pywebview.api.get_state().then((data) => {
    setAppState((previous) => {
      return {
        ...previous,
        ...data,
      };
    });
  });
};

export default useAppState;
