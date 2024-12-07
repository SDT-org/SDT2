import { z } from "zod";
import type { AppState, SetAppState } from "./appState";
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
  };
  raincloud: VisualizationBase & {
    bandwidth: number;
    showMeanline: boolean;
    side: "positive";
    fillColor: ColorString;
    jitter: number;
    markerColor: ColorString;
    markerSize: number;
    pointOpacity: number;
    pointPos: number;
    points: "all" | "outliers" | "suspectedoutliers";
    showAxisLines: boolean;
    showPoints: boolean;
    showZeroLine: boolean;
    plotOrientation: "horizontal" | "vertical";
    editable: boolean;
  };
  violin: VisualizationBase & {
    bandwidth: number;
    boxOpacity: number;
    boxWidth: number;
    boxfillColor: ColorString;
    boxlineColor: ColorString;
    boxlineWidth: number;
    fillColor: ColorString;
    jitter: number;
    markerColor: ColorString;
    markerSize: number;
    pointOpacity: number;
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
  }),
  raincloud: VisualizationBaseSchema.extend({
    bandwidth: z.number(),
    showMeanline: z.boolean(),
    side: z.literal("positive"),
    fillColor: ColorStringSchema,
    jitter: z.number(),
    markerColor: ColorStringSchema,
    markerSize: z.number(),
    pointOpacity: z.number(),
    pointPos: z.number(),
    points: z.enum(["all", "outliers", "suspectedoutliers"]),
    showAxisLines: z.boolean(),
    showPoints: z.boolean(),
    showZeroLine: z.boolean(),
    plotOrientation: z.enum(["horizontal", "vertical"]),
    editable: z.boolean(),
  }),
  violin: VisualizationBaseSchema.extend({
    bandwidth: z.number(),
    boxOpacity: z.number(),
    boxWidth: z.number(),
    boxfillColor: ColorStringSchema,
    boxlineColor: ColorStringSchema,
    boxlineWidth: z.number(),
    fillColor: ColorStringSchema,
    jitter: z.number(),
    markerColor: ColorStringSchema,
    markerSize: z.number(),
    pointOpacity: z.number(),
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
  },
  raincloud: {
    ...visualizationDefaults,
    bandwidth: 8,
    jitter: 0.5,
    markerColor: Colors.Tomato,
    markerSize: 7,
    plotOrientation: "horizontal",
    pointOpacity: 0.5,
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
  },
  violin: {
    ...visualizationDefaults,
    bandwidth: 5,
    boxOpacity: 0.5,
    boxWidth: 0.95,
    boxfillColor: Colors.LightBlue,
    boxlineColor: Colors.Tomato,
    boxlineWidth: 3,
    fillColor: Colors.LightBlue,
    jitter: 0.5,
    markerColor: Colors.Tomato,
    markerSize: 7,
    plotOrientation: "vertical",
    pointOpacity: 0.5,
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
  },
};

export const useDistributionState = (
  appState: AppState,
  setAppState: SetAppState,
) => {
  const updateDistributionState = (values: Partial<DistributionState>) =>
    setAppState((prev) => ({
      ...prev,
      client: {
        ...prev.client,
        distribution: {
          ...prev.client.distribution,
          ...values,
        },
      },
    }));

  const updateVisualization =
    (key: DistributionState["visualization"]) =>
    (
      values: Partial<DistributionState["histogram" | "violin" | "raincloud"]>,
    ) =>
      setAppState((prev) => ({
        ...prev,
        client: {
          ...prev.client,
          distribution: {
            ...prev.client.distribution,
            [key]: {
              ...prev.client.distribution[key],
              ...values,
            },
          },
        },
      }));

  return {
    distributionState: appState.client.distribution,
    updateDistributionState,
    updateHistogram: updateVisualization("histogram"),
    updateRaincloud: updateVisualization("raincloud"),
    updateViolin: updateVisualization("violin"),
  };
};
