import React, { type ErrorInfo } from "react";
import { z } from "zod";
import type { reorderMethods, saveableImageFormats } from "./constants";
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
  type UMAPSettings,
  UMAPSettingsSchema,
} from "./plotTypes";

export const clusterMethodDescriptions = [
  "Unrooted phylogenetic tree for large datasets",
  "Rooted phylogenetic tree ideal for ultra-metric datasets",
];

export type SaveableImageFormat = keyof typeof saveableImageFormats;
export type SaveableImageKey = DocState["dataView"];
export type RasterFormat = Extract<SaveableImageFormat, "png" | "jpeg">;

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
  dataView:
    | "heatmap"
    | "clustermap"
    | "distribution_histogram"
    | "distribution_violin"
    | "umap";
  distribution: DistributionState;
  heatmap: HeatmapSettings;
  clustermap: ClustermapSettings;
  umap: UMAPSettings;
  exportPrefix: string;
  invalid?: {
    reason: string;
  };
  compute_stats?: {
    recommended_cores: number;
    required_memory: number;
    available_memory: number;
  };
  enableClustering: boolean;
  compute_cores?: number;
  cluster_method: keyof typeof reorderMethods;
  analysisMethod: "parasail" | "lzani";
  lzaniScoreType: "ani" | "gani" | "tani";
  overrideParasail: boolean;
  parasail_settings?:
    | {
        scoring_matrix?: string;
        open_penalty?: number;
        extend_penalty?: number;
      }
    | undefined;
  overrideLzani: boolean;
  lzani_settings?:
    | {
        aw?: number;
        am?: number;
        mal?: number;
        msl?: number;
        mrd?: number;
        mqd?: number;
        reg?: number;
        ar?: number;
      }
    | undefined;
  result_metadata?:
    | {
        is_aa?: boolean;
      }
    | undefined;
  exportAlignments: boolean;
  alignmentExportPath: string;
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
  showExportModal: boolean;
  saveFormat: SaveableImageFormat;
  dataExportPath: string;
  lastDataFilePath: string;
  error?: Error | null;
  errorInfo?: ErrorInfo | PromiseRejectionEvent["reason"] | null;
  recentFiles: string[];
  exportStatus: "idle" | "preparing" | "exporting" | "success";
  openExportFolder: boolean;
};

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
    "umap",
  ]),
  distribution: DistributionStateSchema,
  heatmap: HeatmapSettingsSchema,
  clustermap: ClustermapSettingsSchema,
  umap: UMAPSettingsSchema,
});

const clientDocState = {
  exportPrefix: "",
  openExportFolder: false,
};

export const initialDocState: DocState = {
  id: "",
  view: "runner",
  filename: "",
  filetype: "",
  filemtime: null,
  basename: "",
  modified: false,
  parsed: false,
  progress: 0,
  sequences_count: 0,
  stage: "Preprocessing",
  pair_progress: 0,
  pair_count: 0,
  dataView: "heatmap",
  distribution: initialDistributionState,
  enableClustering: true,
  analysisMethod: "parasail",
  lzaniScoreType: "ani",
  cluster_method: "average",
  overrideParasail: false,
  overrideLzani: false,
  exportAlignments: false,
  alignmentExportPath: "",
  heatmap: {
    colorScaleKey: "Portland",
    reverse: false,
    vmax: 100,
    vmin: 65,
    cellspace: 1,
    annotation: false,
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
    axlabel_fontsize: 12,
    axlabel_yrotation: 0,
    cutoff_1: 95,
    cutoff_2: 75,
    hideValuesBelow: 0,
    hideValuesBelowEnabled: false,
  },
  clustermap: {
    threshold: 70,
    method: "average",
    annotation: true,
    titleFont: "Sans Serif",
    showTitles: true,
    title: "",
    xtitle: "",
    ytitle: "",
    axis_labels: false,
    axlabel_xrotation: 0,
    axlabel_fontsize: 12,
    axlabel_yrotation: 0,
    cellspace: 1,
    showLegend: true,
    showClusterCounts: false,
  },
  umap: {
    n_neighbors: 15,
    min_dist: 0.1,
    pointSize: 4,
    opacity: 0.8,
  },
  ...clientDocState,
};

export const initialAppState: AppState = {
  activeDocumentId: "",
  debug: false,
  dataExportPath: "",
  lastDataFilePath: "",
  saveFormat: "svg",
  showExportModal: false,
  documents: [],
  platform: {
    cores: 1,
    memory: 1,
    platform: "unknown",
  },
  recentFiles: [],
  exportStatus: "idle",
  openExportFolder: false,
};

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
