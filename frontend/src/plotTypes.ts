import { z } from "zod";
import type { colorScales } from "./colorScales";

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
});

export type HeatmapData = Array<Array<number | null>>;
export type MetaData = GetDataResponse["metadata"];

export type GetDataResponse = {
  data: HeatmapData & string[][];
  metadata: {
    minVal: number;
    maxVal: number;
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
  method: "single" | "complete" | "average" | "weighted";
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
}

export const ClustermapSettingsSchema = z.object({
  threshold_: z.number(),
  method: z.enum(["single", "complete", "average", "weighted"]),
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
});
