import React, { type ErrorInfo } from "react";
import {
  type DistributionState,
  initialDistributionState,
} from "./distributionState";
import type messages from "./messages";
import type { HeatmapSettings } from "./plotTypes";

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
    dataView: "heatmap" | "distribution";
    enableClustering: boolean;
    enableOutputAlignments: boolean;
    cluster_method: (typeof clusterMethods)[number];
    compute_cores: number;
    error?: Error | null;
    errorInfo?: ErrorInfo | PromiseRejectionEvent["reason"] | null;
    saveFormat: SaveableImageFormat;
    showExportModal: boolean;
    alignmentExportPath: string;
    dataExportPath: string;
    lastDataFilePath: string;
    distribution: DistributionState;
    heatmap: HeatmapSettings;
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
  stage: "Preprocessing",
  pair_progress: 0,
  pair_count: 0,
  client: {
    dataView: "heatmap",
    dataExportPath: "",
    alignmentExportPath: "",
    lastDataFilePath: "",
    enableClustering: true,
    enableOutputAlignments: false,
    cluster_method: "Neighbor-Joining",
    saveFormat: "svg",
    showExportModal: false,
    compute_cores: 1,
    distribution: initialDistributionState,
    heatmap: {
      colorscale: "Portland",
      reverse: false,
      vmax: 100,
      vmin: 65,
      cellspace: 1,
      annotation: false,
      annotation_font_size: 10,
      annotation_rounding: 0,
      annotation_alpha: "0",
      color: "white",
      showscale: true,
      cbar_shrink: 1,
      cbar_aspect: 25,
      cbar_pad: 10,
      axis_labels: false,
      axlabel_xrotation: 270,
      axlabel_xfontsize: 12,
      axlabel_yrotation: 360,
      axlabel_yfontsize: 12,
      cutoff_1: 95,
      cutoff_2: 75,
    },
  },
  platform: {
    cores: 1,
    memory: 1,
    platform: "unknown",
  },
};

export const clientStateKey = "app-client-state";

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
