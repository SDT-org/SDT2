import { z } from "zod";
import type { colorScales } from "./colorScales";

export type Colorscale = keyof typeof colorScales;

export interface HeatmapSettings {
  colorscale: Colorscale | "Discrete";
  reverse: boolean;
  vmax: number;
  vmin: number;
  cellspace: number;
  annotation: boolean;
  annotation_font_size: number;
  annotation_rounding: 0 | 1 | 2;
  annotation_alpha: string;
  showscale: boolean;
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
  colorscale: z.enum([
    "Greys",
    "YlGnBu",
    "Greens",
    "YlOrRd",
    "Bluered",
    "RdBu",
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

export type HeatmapData = string[][];

export type GetDataResponse = {
  data: string[][];
  metadata: {
    minVal: number;
    maxVal: number;
  };
  gc_stats: number[];
  length_stats: number[];
  identity_scores: [string, string, number][];
};

export type DistributionData = Omit<
  GetDataResponse,
  "data" | "identity_scores"
> & {
  raw_mat: number[];
  identity_combos: [string, string][];
};
