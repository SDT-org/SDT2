import { z } from "zod";
import type { DocState, SetDocState } from "./appState";
import { type ColorString, ColorStringSchema, Colors } from "./colors";

export type Visualization = "histogram" | "violin" | "raincloud";
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
  visualization: Visualization;
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
    subtitle: string;
    xtitle: string;
    ytitle: string;
    plotOrientation: "horizontal" | "vertical";
    titleFont: "Monospace" | "Sans Serif";
  };
  raincloud: VisualizationBase & {
    bandwidth: number;
    showMeanline: boolean;
    side: "positive";
    fillColor: ColorString;
    jitter: number;
    markerColor: ColorString;
    markerSize: number;
    pointPos: number;
    points: "all" | "outliers" | "suspectedoutliers";
    showAxisLines: boolean;
    showPoints: boolean;
    showZeroLine: boolean;
    plotOrientation: "horizontal" | "vertical";
    editable: boolean;
    dticks: number;
    title: string;
    subtitle: string;
    xtitle: string;
    ytitle: string;
    titleFont: "Monospace" | "Sans Serif";
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
    subtitle: string;
    xtitle: string;
    ytitle: string;
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
  visualization: z.enum(["histogram", "violin", "raincloud"]),
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
    subtitle: z.string(),
    xtitle: z.string(),
    ytitle: z.string(),
    plotOrientation: z.enum(["horizontal", "vertical"]),
    titleFont: z.enum(["Monospace", "Sans Serif"]),
  }),
  raincloud: VisualizationBaseSchema.extend({
    bandwidth: z.number(),
    showMeanline: z.boolean(),
    side: z.literal("positive"),
    fillColor: ColorStringSchema,
    jitter: z.number(),
    markerColor: ColorStringSchema,
    markerSize: z.number(),
    pointPos: z.number(),
    points: z.enum(["all", "outliers", "suspectedoutliers"]),
    showAxisLines: z.boolean(),
    showPoints: z.boolean(),
    showZeroLine: z.boolean(),
    plotOrientation: z.enum(["horizontal", "vertical"]),
    editable: z.boolean(),
    dticks: z.number(),
    title: z.string(),
    subtitle: z.string(),
    xtitle: z.string(),
    ytitle: z.string(),
    titleFont: z.enum(["Monospace", "Sans Serif"]),
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
    subtitle: z.string(),
    xtitle: z.string(),
    ytitle: z.string(),
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
  visualization: "histogram",
  dataSet: "scores",
  histogram: {
    ...visualizationDefaults,
    binColor: Colors.LightBlue,
    binSize: 1,
    histOutlineWidth: 1,
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
    subtitle: "Histogram",
    title: "Histogram",
    xtitle: "Percent Identity",
    ytitle: "Frequency",
    plotOrientation: "horizontal",
    titleFont: "Sans Serif",
  },
  raincloud: {
    ...visualizationDefaults,
    bandwidth: 8,
    jitter: 0.5,
    markerColor: Colors.Tomato,
    markerSize: 7,
    plotOrientation: "horizontal",
    pointPos: -1.5,
    points: "all",
    showAxisLines: true,
    showPoints: true,
    showZeroLine: false,
    fillColor: Colors.LightBlue,
    editable: false,
    side: "positive",
    showMeanline: true,
    makeEditable: true,
    dticks: 5,
    showTitles: true,
    subtitle: "Raincloud Plot",
    title: "Raincloud Plot",
    xtitle: "Percent Identity",
    ytitle: "Genome",
    titleFont: "Sans Serif",
  },
  violin: {
    ...visualizationDefaults,
    bandwidth: 5,
    boxWidth: 0.5,
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
    subtitle: "Violin Plot",
    xtitle: "",
    ytitle: "",
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
    (key: DistributionState["visualization"]) =>
    (
      values: Partial<DistributionState["histogram" | "violin" | "raincloud"]>,
    ) =>
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
    updateRaincloud: updateVisualization("raincloud"),
    updateViolin: updateVisualization("violin"),
  };
};
