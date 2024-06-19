import React, { ErrorInfo } from "react";

export const performanceProfiles = ["balanced", "best", "low"] as const;
export const alignmentTypes = ["global", "local"] as const;
export const clusterMethods = ["Neighbor-Joining", "UPGMA", "None"] as const;
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
  client: {
    dataView: "heatmap" | "plot";
    cluster_method: (typeof clusterMethods)[keyof typeof clusterMethods];
    performance_profile: (typeof performanceProfiles)[keyof typeof performanceProfiles];
    alignment_type: (typeof alignmentTypes)[keyof typeof alignmentTypes];
    error?: Error;
    errorInfo?: ErrorInfo;
    saveFormat: "png" | "jpeg" | "svg";
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
  client: {
    dataView: "heatmap",
    cluster_method: "Neighbor-Joining",
    performance_profile: "best",
    alignment_type: "global",
    saveFormat: "svg",
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
    return;
  }

  return window.pywebview.api.get_state().then((data) => {
    setAppState((previous) => {
      return {
        ...previous,
        ...data,
      };
    });
  });

export default useAppState;
