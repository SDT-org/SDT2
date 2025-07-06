import { z } from "zod";
import type { DocState, SetDocState } from "./appState";
import { type ColorString, ColorStringSchema, Colors } from "./colors";

export type Visualization = "distribution_histogram" | "distribution_violin";
type DataSet = number[];
export type DataSets = {
  scores: DataSet;
  gc: DataSet;
  length: DataSet;
};

type VisualizationBase = {
  plotTitle: string;
  lineColor: ColorString;
  lineWidth: number;
  showGrid: boolean;
  showTickLabels: boolean;
  showAxisLabels: boolean;
  showAxisLines: boolean;
  makeEditable: boolean;
  showMeanline: boolean;
  showTitles: boolean;
};

export type DistributionState = {
  dataSet: keyof DataSets;
  histogram: VisualizationBase & {
    binColor: ColorString;
    binSize: number;
    histOutlineWidth: number;
    histnorm: "probability" | "percent";
    showHistogram: boolean;
    showLine: boolean;
    histlineColor: ColorString;
    dtickx: number;
    dticky: number;
    title: string;
    plotOrientation: "horizontal" | "vertical";
    titleFont: "Monospace" | "Sans Serif";
    barGap: number;
  };
  violin: VisualizationBase & {
    bandwidth: number;
    boxWidth: number;
    boxfillColor: ColorString;
    boxlineColor: ColorString;
    boxlineWidth: number;
    fillColor: ColorString;
    jitter: number;
    markerColor: ColorString;
    markerSize: number;
    pointOrientation: "Violin" | "Box";
    pointPos: number;
    points: "all" | "outliers" | "suspectedoutliers";
    showAxisLines: boolean;
    showBox: boolean;
    showMeanline: boolean;
    showPoints: boolean;
    showViolin: boolean;
    showZeroLine: boolean;
    whiskerWidth: number;
    plotOrientation: "horizontal" | "vertical";
    title: string;
    titleFont: "Monospace" | "Sans Serif";
  };
};

const VisualizationBaseSchema = z.object({
  plotTitle: z.string(),
  lineColor: ColorStringSchema,
  lineWidth: z.number(),
  showGrid: z.boolean(),
  showTickLabels: z.boolean(),
  showAxisLabels: z.boolean(),
  showAxisLines: z.boolean(),
  makeEditable: z.boolean(),
  showMeanline: z.boolean(),
  showTitles: z.boolean(),
});

export const DistributionStateSchema = z.object({
  dataSet: z.enum(["scores", "gc", "length"]),
  histogram: VisualizationBaseSchema.extend({
    binColor: ColorStringSchema,
    binSize: z.number(),
    histOutlineWidth: z.number(),
    histnorm: z.enum(["probability", "percent"]),
    showHistogram: z.boolean(),
    showLine: z.boolean(),
    histlineColor: ColorStringSchema,
    dtickx: z.number(),
    dticky: z.number(),
    title: z.string(),
    plotOrientation: z.enum(["horizontal", "vertical"]),
    titleFont: z.enum(["Monospace", "Sans Serif"]),
    barGap: z.number(),
  }),
  violin: VisualizationBaseSchema.extend({
    bandwidth: z.number(),
    boxWidth: z.number(),
    boxfillColor: ColorStringSchema,
    boxlineColor: ColorStringSchema,
    boxlineWidth: z.number(),
    fillColor: ColorStringSchema,
    jitter: z.number(),
    markerColor: ColorStringSchema,
    markerSize: z.number(),
    pointOrientation: z.enum(["Violin", "Box"]),
    pointPos: z.number(),
    points: z.enum(["all", "outliers", "suspectedoutliers"]),
    showAxisLines: z.boolean(),
    showBox: z.boolean(),
    showMeanline: z.boolean(),
    showPoints: z.boolean(),
    showViolin: z.boolean(),
    showZeroLine: z.boolean(),
    whiskerWidth: z.number(),
    plotOrientation: z.enum(["horizontal", "vertical"]),
    title: z.string(),
    titleFont: z.enum(["Monospace", "Sans Serif"]),
  }),
});

const visualizationDefaults = {
  plotTitle: "Distribution of Percent Identities",
  lineColor: Colors.Tomato,
  lineWidth: 3,
  showAxisLabels: true,
  showGrid: true,
  showTickLabels: true,
};

export const initialDistributionState: DistributionState = {
  dataSet: "scores",
  histogram: {
    ...visualizationDefaults,
    binColor: Colors.LightBlue,
    binSize: 1,
    histOutlineWidth: 0,
    histlineColor: Colors.Tomato,
    histnorm: "probability",
    showHistogram: true,
    showLine: true,
    makeEditable: true,
    showAxisLines: true,
    showMeanline: true,
    dtickx: 5,
    dticky: 1,
    showTitles: true,
    title: "Histogram",
    plotOrientation: "horizontal",
    titleFont: "Sans Serif",
    barGap: 0.1,
  },
  violin: {
    ...visualizationDefaults,
    bandwidth: 5,
    boxWidth: 0.05,
    boxfillColor: Colors.LightBlue,
    boxlineColor: Colors.Tomato,
    boxlineWidth: 3,
    fillColor: Colors.LightBlue,
    jitter: 0.5,
    markerColor: Colors.Tomato,
    markerSize: 7,
    plotOrientation: "vertical",
    pointOrientation: "Violin",
    pointPos: 0,
    points: "all",
    showAxisLines: true,
    showBox: true,
    showMeanline: true,
    makeEditable: true,
    showPoints: true,
    showViolin: true,
    showZeroLine: false,
    whiskerWidth: 0.2,
    showTitles: true,
    showTickLabels: true,
    title: "Violin Plot",
    titleFont: "Sans Serif",
  },
};

export const useDistributionState = (
  docState: DocState,
  setDocState: SetDocState,
) => {
  const updateDistributionState = (values: Partial<DistributionState>) =>
    setDocState((prev) => ({
      ...prev,
      distribution: {
        ...prev.distribution,
        ...values,
      },
    }));

  const updateVisualization =
    (key: "histogram" | "violin") =>
    (values: Partial<DistributionState["histogram" | "violin"]>) =>
      setDocState((prev) => ({
        ...prev,
        distribution: {
          ...prev.distribution,
          [key]: {
            ...prev.distribution[key],
            ...values,
          },
        },
      }));

  return {
    distributionState: docState.distribution,
    updateDistributionState,
    updateHistogram: updateVisualization("histogram"),
    updateViolin: updateVisualization("violin"),
  };
};
