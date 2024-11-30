import React, { type ErrorInfo } from "react";
import type messages from "./messages";

export const clusterMethods = ["Neighbor-Joining", "UPGMA"] as const;

export const clusterMethodDescriptions = [
  "Unrooted phylogenetic tree for large datasets",
  "Rooted phylogenetic tree ideal for ultra-metric datasets",
];

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
    enableClustering: boolean;
    enableOutputAlignments: boolean;
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
    enableClustering: true,
    enableOutputAlignments: false,
    cluster_method: "Neighbor-Joining",
    saveFormat: "svg",
    showExportModal: false,
    compute_cores: 1,
  },
  platform: {
    cores: 1,
    memory: 1,
    platform: "unknown",
  },
};

export type SetAppState = React.Dispatch<React.SetStateAction<AppState>>;
export type SyncStateEvent = CustomEvent<{ state: AppState }>;

export const AppStateContext = React.createContext<{
  appState: AppState;
  setAppState: SetAppState;
}>({
  appState: initialAppState,
  setAppState: () => {
    throw new Error("setAppState must be initialized");
  },
});

export const useAppState = () => {
  const context = React.useContext(AppStateContext);
  if (context === undefined) {
    throw new Error("useAppState must be used within an AppStateProvider");
  }
  return context;
};

export default useAppState;
