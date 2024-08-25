import React, { ErrorInfo } from "react";
import messages from "./messages";

export const performanceProfiles = {
  best: "Best",
  high: "High Performance",
  balanced: "Balanced",
  low: "Energy Saver",
} as const;
export const clusterMethods = ["Neighbor-Joining", "UPGMA", "None"] as const;

export type SaveableImageFormat = "png" | "jpeg" | "svg";

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
  performance_profiles: { [_: string]: number };
  stage: string;
  pair_progress: number;
  pair_count: number;
  estimated_time?: number;
  validation_error_id?: keyof typeof messages;
  client: {
    dataView: "heatmap" | "plot";
    cluster_method: (typeof clusterMethods)[number];
    performance_profile: keyof typeof performanceProfiles;
    error?: Error;
    errorInfo?: ErrorInfo;
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
  // These are just to make frontend easier to test, they get overwritten during the initial syncAppState
  performance_profiles: { best: 4, high: 3, balanced: 2, low: 1 },
  stage: "Preprocessing",
  pair_progress: 0,
  pair_count: 0,
  client: {
    dataView: "heatmap",
    cluster_method: "Neighbor-Joining",
    performance_profile: "best",
    saveFormat: "svg",
    showExportModal: false,
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
