import React, { type ErrorInfo } from "react";
import { z } from "zod";
import {
  type DistributionState,
  DistributionStateSchema,
  initialDistributionState,
} from "./distributionState";
import type messages from "./messages";
import {
  type ClustermapSettings,
  ClustermapSettingsSchema,
  type HeatmapSettings,
  HeatmapSettingsSchema,
} from "./plotTypes";
export const clusterMethods = ["Neighbor-Joining", "UPGMA"] as const;

export const clusterMethodDescriptions = [
  "Unrooted phylogenetic tree for large datasets",
  "Rooted phylogenetic tree ideal for ultra-metric datasets",
];

export const saveableImageFormats = {
  svg: "SVG",
  png: "PNG",
  jpeg: "JPEG",
};
export type SaveableImageFormat = keyof typeof saveableImageFormats;
const saveableImageFormatKeys = ["svg", "png", "jpeg"] as const;

export type DocState = {
  id: string;
  view: "runner" | "loader" | "viewer";
  filename: string;
  filetype: string;
  filemtime: number | null;
  basename: string;
  modified: boolean;
  parsed: boolean;
  progress: number;
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
  dataView:
    | "heatmap"
    | "clustermap"
    | "distribution_histogram"
    | "distribution_violin"
    | "distribution_raincloud";
  distribution: DistributionState;
  heatmap: HeatmapSettings;
  clustermap: ClustermapSettings;
};

export type AppState = {
  activeDocumentId: DocState["id"];
  active_run_document_id?: DocState["id"];
  documents: DocState[];
  // TODO: make platform info an api, doesn't need to be state at all
  platform: {
    cores: number;
    memory: number;
    platform: string;
  };
  config?: {
    appVersion: string;
    userPath: string;
  };
  debug: boolean;
  enableClustering: boolean;
  enableOutputAlignments: boolean;
  cluster_method: (typeof clusterMethods)[number];
  compute_cores: number;
  showExportModal: boolean;
  saveFormat: SaveableImageFormat;
  alignmentExportPath: string;
  dataExportPath: string;
  lastDataFilePath: string;
  error?: Error | null;
  errorInfo?: ErrorInfo | PromiseRejectionEvent["reason"] | null;
};

export const clientStateSchema = z.object({
  enableClustering: z.boolean(),
  enableOutputAlignments: z.boolean(),
  cluster_method: z.enum(clusterMethods),
  compute_cores: z.number(),
  error: z.instanceof(Error).nullable(),
  errorInfo: z.null(),
  saveFormat: z.enum(saveableImageFormatKeys),
  showExportModal: z.boolean(),
  alignmentExportPath: z.string(),
  dataExportPath: z.string(),
  lastDataFilePath: z.string(),
});

export const docStateSchema = z.object({
  id: z.string(),
  view: z.enum(["runner", "loader", "viewer"]),
  filename: z.string(),
  filetype: z.string(),
  filemtime: z.number().nullable(),
  basename: z.string(),
  modified: z.boolean(),
  progress: z.number(),
  sequences_count: z.number(),
  pair_progress: z.number(),
  pair_count: z.number(),
  stage: z.enum([
    "",
    "Preparing",
    "Preprocessing",
    "Analyzing",
    "Postprocessing",
    "Finalizing",
    "Processed",
  ]),
  dataView: z.enum([
    "heatmap",
    "clustermap",
    "distribution_histogram",
    "distribution_violin",
    "distribution_raincloud",
  ]),
  distribution: DistributionStateSchema,
  heatmap: HeatmapSettingsSchema,
  clustermap: ClustermapSettingsSchema,
});

export const initialDocState: DocState = {
  id: "",
  view: "runner",
  filename: "",
  filetype: "",
  filemtime: null,
  basename: "",
  modified: false,
  parsed: false,
  // TODO: wrap in runState
  progress: 0,
  sequences_count: 0,
  stage: "Preprocessing",
  pair_progress: 0,
  pair_count: 0,
  dataView: "heatmap",
  distribution: initialDistributionState,
  heatmap: {
    colorScaleKey: "Portland",
    reverse: false,
    vmax: 100,
    vmin: 65,
    cellspace: 1,
    annotation: false,
    annotation_font_size: 10,
    annotation_rounding: 0,
    showscale: true,
    titleFont: "Sans Serif",
    showTitles: false,
    title: "",
    xtitle: "",
    ytitle: "",
    cbar_shrink: 5,
    cbar_aspect: 2.5,
    cbar_pad: 10,
    axis_labels: false,
    axlabel_xrotation: 0,
    axlabel_xfontsize: 12,
    axlabel_yrotation: 0,
    axlabel_yfontsize: 12,
    cutoff_1: 95,
    cutoff_2: 75,
  },
  clustermap: {
    threshold_one: 85,
    threshold_two: 0,
    annotation: false,
    annotation_font_size: 10,
    annotation_rounding: 0,
    titleFont: "Sans Serif",
    showTitles: false,
    title: "",
    xtitle: "",
    ytitle: "",
    axis_labels: false,
    axlabel_xrotation: 0,
    axlabel_xfontsize: 12,
    axlabel_yrotation: 0,
    axlabel_yfontsize: 12,
    cellspace: 1,
  },
};

export const initialAppState: AppState = {
  activeDocumentId: "",
  debug: false,
  dataExportPath: "",
  alignmentExportPath: "",
  lastDataFilePath: "",
  enableClustering: true,
  enableOutputAlignments: false,
  cluster_method: "Neighbor-Joining",
  saveFormat: "svg",
  showExportModal: false,
  compute_cores: 1,
  documents: [],
  platform: {
    cores: 1,
    memory: 1,
    platform: "unknown",
  },
};

export const clientStateKey = "app-client-state";

export type SetAppState = React.Dispatch<React.SetStateAction<AppState>>;
export type SetDocState = (
  nextDoc: (prevDoc: DocState) => DocState,
  markModified?: boolean,
) => void;
export type UpdateDocState = (
  newValues: Partial<DocState>,
  markModified?: boolean,
) => void;
export type SyncStateEvent = CustomEvent<{ state: AppState }>;
export type SyncProgressEvent = CustomEvent<
  Pick<DocState, "id" | "progress" | "pair_progress" | "estimated_time">
>;

export const findDoc = (id: string, documents: AppState["documents"]) =>
  documents.find((d) => d.id === id);

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

export type HeatmapRefType = HTMLCanvasElement | SVGSVGElement | null;
export const HeatmapRefContext =
  React.createContext<React.MutableRefObject<HeatmapRefType> | null>(null);
