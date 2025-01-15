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
  annotation_font_size: number;
  annotation_rounding: 0 | 1 | 2;
  annotation_alpha: string;
  showscale: boolean;
  titleFont: "Sans Serif" | "Monospace";
  showTitles: boolean;
  title: string;
  subtitle: string;
  xtitle: string;
  ytitle: string;
  cbar_shrink: number;
  cbar_pad: number;
  cbar_aspect: number;
  axis_labels: boolean;
  axlabel_xrotation: number;
  axlabel_xfontsize: number;
  axlabel_yrotation: number;
  axlabel_yfontsize: number;
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
  annotation_font_size: z.number(),
  annotation_rounding: z.union([z.literal(0), z.literal(1), z.literal(2)]),
  annotation_alpha: z.string(),
  showscale: z.boolean(),
  showTitles: z.boolean(),
  titleFont: z.enum(["Sans Serif", "Monospace"]),
  title: z.string(),
  subtitle: z.string(),
  xtitle: z.string(),
  ytitle: z.string(),
  cbar_shrink: z.number(),
  cbar_pad: z.number(),
  cbar_aspect: z.number(),
  axis_labels: z.boolean(),
  axlabel_xrotation: z.number(),
  axlabel_xfontsize: z.number(),
  axlabel_yrotation: z.number(),
  axlabel_yfontsize: z.number(),
  cutoff_1: z.number(),
  cutoff_2: z.number(),
});

export type HeatmapData = GetDataResponse["data"];
export type MetaData = GetDataResponse["metadata"];

export type GetDataResponse = {
  data: string[][];
  metadata: {
    minVal: number;
    maxVal: number;
  };
  identity_scores: [string, string, number][];
  stat_ids: string[];
  full_stats: [string, number, number][];
};

// TODO: make this file about heatmap
export type DistributionData = Omit<
  GetDataResponse,
  "data" | "identity_scores" | "metadata" | "stat_ids"
> & {
  raw_mat: number[];
  identity_combos: [string, string][];
  gc_stats: number[];
  length_stats: number[];
};
