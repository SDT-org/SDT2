import { z } from "zod";
import type { colorScales } from "./colorScales";
import { reorderMethods } from "./constants";

export type ColorScaleKey = keyof typeof colorScales;

export interface HeatmapSettings {
  colorScaleKey: ColorScaleKey | "Discrete";
  reverse: boolean;
  vmax: number;
  vmin: number;
  cellspace: number;
  annotation: boolean;
  annotation_rounding: 0 | 1 | 2;
  showscale: boolean;
  titleFont: "Sans Serif" | "Monospace";
  showTitles: boolean;
  title: string;
  xtitle: string;
  ytitle: string;
  cbar_shrink: number;
  cbar_pad: number;
  cbar_aspect: number;
  axis_labels: boolean;
  axlabel_xrotation: number;
  axlabel_fontsize: number;
  axlabel_yrotation: number;
  cutoff_1: number;
  cutoff_2: number;
  hideValuesBelow: number;
  hideValuesBelowEnabled: boolean;
}

export const HeatmapSettingsSchema = z.object({
  colorScaleKey: z.enum([
    "Greys",
    "Greens",
    "Yellow_Green_Blue",
    "Yellow_Orange_Red",
    "Blue_Red",
    "Red_Blue",
    "Reds",
    "Blues",
    "Picnic",
    "Rainbow",
    "Portland",
    "Jet",
    "Hot",
    "Blackbody",
    "Earth",
    "Electric",
    "Viridis",
    "Cividis",
    "Discrete",
  ]),
  reverse: z.boolean(),
  vmax: z.number(),
  vmin: z.number(),
  cellspace: z.number(),
  annotation: z.boolean(),
  annotation_rounding: z.union([z.literal(0), z.literal(1), z.literal(2)]),
  showscale: z.boolean(),
  showTitles: z.boolean(),
  titleFont: z.enum(["Sans Serif", "Monospace"]),
  title: z.string(),
  xtitle: z.string(),
  ytitle: z.string(),
  cbar_shrink: z.number(),
  cbar_pad: z.number(),
  cbar_aspect: z.number(),
  axis_labels: z.boolean(),
  axlabel_xrotation: z.number(),
  axlabel_fontsize: z.number(),
  axlabel_yrotation: z.number(),
  cutoff_1: z.number(),
  cutoff_2: z.number(),
  hideValuesBelow: z.number(),
  hideValuesBelowEnabled: z.boolean(),
});

export type HeatmapData = Array<Array<number | null>>;
export type MetaData = GetDataResponse["metadata"];

export type GetDataResponse = {
  data: HeatmapData & string[][];
  metadata: {
    minVal: number;
    maxVal: number;
    run?: {
      analysis_method: "parasail" | "lzani";
      cluster_method: (typeof reorderMethodKeys)[number];
      lzani: {
        score_type: string;
      };
    };
    unaligned_count?: number;
    distribution_stats?: {
      mean: number;
      median: number;
      std: number;
      min: number;
      max: number;
      q1: number;
      q3: number;
      count: number;
    };
    gc_stats?: {
      mean: number;
      median: number;
      std: number;
      min: number;
      max: number;
      q1: number;
      q3: number;
      count: number;
    };
    length_stats?: {
      mean: number;
      median: number;
      std: number;
      min: number;
      max: number;
      q1: number;
      q3: number;
      count: number;
    };
  };
  ids: string[];
  identity_scores: [number, number, number][];
  stat_ids: string[];
  full_stats: [string, number, number][];
};

// TODO: make this file about heatmap
export type DistributionData = Omit<
  GetDataResponse,
  "data" | "identity_scores" | "metadata" | "stat_ids"
> & {
  raw_mat: number[];
  ids: string[];
  identity_combos: [number, number][];
  gc_stats: number[];
  length_stats: number[];
};

export interface ClustermapSettings {
  threshold: number;
  method: keyof typeof reorderMethods;
  annotation: boolean;
  titleFont: "Sans Serif" | "Monospace";
  showTitles: boolean;
  title: string;
  xtitle: string;
  ytitle: string;
  axis_labels: boolean;
  axlabel_xrotation: number;
  axlabel_fontsize: number;
  axlabel_yrotation: number;
  cellspace: number;
  showLegend: boolean;
  showClusterCounts: boolean;
}

const reorderMethodKeys = Object.keys(
  reorderMethods,
) as (keyof typeof reorderMethods)[];

export const ClustermapSettingsSchema = z.object({
  threshold: z.number(),
  method: z.any().refine((v) => reorderMethodKeys.includes(v)), // TODO:fix this
  annotation: z.boolean(),
  showTitles: z.boolean(),
  titleFont: z.enum(["Sans Serif", "Monospace"]),
  title: z.string(),
  xtitle: z.string(),
  ytitle: z.string(),
  axis_labels: z.boolean(),
  axlabel_xrotation: z.number(),
  axlabel_fontsize: z.number(),
  axlabel_yrotation: z.number(),
  cellspace: z.number(),
  showLegend: z.boolean(),
  showClusterCounts: z.boolean(),
});

export type ClusterDataItem = {
  id: string;
  cluster: number;
  original_cluster?: number;
};

export type ClusterStats = {
  total_clusters: number;
  largest_cluster: number;
  singleton_clusters: number;
};

export type GetClustermapDataResponse = {
  matrix: HeatmapData;
  tickText: string[];
  clusterData: ClusterDataItem[];
  cluster_stats: ClusterStats;
};

export interface UMAPSettings {
  n_neighbors: number;
  min_dist: number;
  pointSize: number;
  opacity: number;
  minClusterSize: number;
  clusterEpsilon: number;
  showClusterBoundaries: boolean;
  colorByCluster: boolean;
  colorBy: "cluster" | "metadata";
  selectedMetadataColumn?: string;
  uploadedMetadata?: {
    columns: string[];
    columnTypes: Record<string, string>;
    matchStats: {
      totalMetadataIds: number;
      totalSequenceIds: number;
      exactMatches: number;
      versionMatches: number;
      unmatched: number;
      matchPercentage: number;
    };
  };
}

export const UMAPSettingsSchema = z.object({
  n_neighbors: z.number().min(2).max(200),
  min_dist: z.number().min(0).max(1),
  pointSize: z.number().min(1).max(20),
  opacity: z.number().min(0.1).max(1),
  minClusterSize: z.number().min(2).max(50),
  clusterEpsilon: z.number().min(0).max(50),
  showClusterBoundaries: z.boolean(),
  colorByCluster: z.boolean(),
  colorBy: z.enum(["cluster", "metadata"]),
  selectedMetadataColumn: z.string().optional(),
  uploadedMetadata: z
    .object({
      columns: z.array(z.string()),
      columnTypes: z.record(z.string()),
      matchStats: z.object({
        totalMetadataIds: z.number(),
        totalSequenceIds: z.number(),
        exactMatches: z.number(),
        versionMatches: z.number(),
        unmatched: z.number(),
        matchPercentage: z.number(),
      }),
    })
    .optional(),
});

export type UMAPPoint = {
  id: string;
  x: number;
  y: number;
  cluster?: number;
  clusters?: {
    [method: string]: number;
  };
};

export type UMAPData = {
  embedding: UMAPPoint[];
  bounds: {
    x: [number, number];
    y: [number, number];
  };
  clusterMethods?: {
    [method: string]: { [id: string]: number };
  };
  threshold?: number;
};

export type GetUMAPDataResponse = {
  data: UMAPData;
  metadata: {
    n_neighbors: number;
    min_dist: number;
    sequences_count: number;
    cluster_method?: string;
    cluster_methods?: string[];
    threshold?: number;
    min_cluster_size?: number;
    cluster_epsilon?: number;
  };
};
