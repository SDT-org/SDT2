import tinycolor from "tinycolor2";
import type { ColorScaleArray } from "./colorScales";
import { createD3ColorScale, distinctColor } from "./colors";
import type {
  ClusterDataItem,
  HeatmapData,
  HeatmapSettings,
} from "./plotTypes";

export const getPlotMetrics = (width: number, height: number) => {
  // WIP
  return {
    width: width,
    height: height,
  };
};

export const getCellMetrics = (
  cellSize: number,
  cellSpace: number,
  characterCount: number,
) => {
  if (characterCount <= 0) throw new Error("characterCount must be > 0");

  const CHARACTER_WIDTH = 0.6; // Approximate width for Roboto Mono @ 10px
  const USABLE_SPACE_RATIO = 0.8;
  const MIN_FONT_SIZE = 0.25;
  const MAX_FONT_SIZE = 20;

  const scaledCellSpace = cellSpace * (cellSize / (cellSize + 20));
  const spacedCellSize = Math.max(1, cellSize - scaledCellSpace);
  const scaledFontSize = Math.min(
    MAX_FONT_SIZE,
    Math.max(
      MIN_FONT_SIZE,
      (spacedCellSize * USABLE_SPACE_RATIO) /
        (characterCount * CHARACTER_WIDTH),
    ),
  );

  return {
    cellSize: spacedCellSize,
    cellSpace: scaledCellSpace,
    cellOffset: scaledCellSpace / 2,
    fontSize: scaledFontSize,
    textOffset: spacedCellSize / 2,
  };
};

export const formatHeatmapData = (
  data: HeatmapData,
  settings: Pick<
    HeatmapSettings,
    "colorScaleKey" | "vmax" | "vmin" | "annotation_rounding"
  >,
  colorScale: ColorScaleArray,
) => {
  const colorFn = createD3ColorScale(
    colorScale,
    settings.colorScaleKey === "Discrete",
    settings.vmax,
    settings.vmin,
  );

  return data.flatMap((row, y) =>
    row.filter(Number).map((value, x) => {
      const backgroundColor = colorFn(Number(value));
      const foregroundColor = tinycolor(backgroundColor).isLight()
        ? "#000"
        : "#fff";
      const roundedValue = (value as number).toFixed(
        settings.annotation_rounding,
      );

      return {
        x,
        y,
        value: value as number,
        displayValue: value === 100 ? "100" : roundedValue.toString(),
        backgroundColor,
        foregroundColor,
      };
    }),
  );
};

export const formatClustermapData = (
  data: HeatmapData,
  tickText: string[],
  clusterData?: ClusterDataItem[],
) => {
  const clusterMap = new Map();
  if (clusterData) {
    for (const item of clusterData) {
      clusterMap.set(item.id, item.cluster);
    }
  }

  const colorCache = new Map();

  const result = data.flatMap((row, y) =>
    row.filter(Number).map((value, x) => {
      const clusterX = clusterMap.get(tickText[x]);
      const clusterY = clusterMap.get(tickText[y]);

      const clusterMatch =
        clusterX !== undefined &&
        clusterY !== undefined &&
        clusterX === clusterY;

      const clusterGroup = clusterMatch ? clusterX : null;

      let backgroundColor = "rgb(245, 245, 245)";
      if (clusterGroup) {
        if (!colorCache.has(clusterGroup)) {
          colorCache.set(clusterGroup, distinctColor(clusterGroup));
        }
        backgroundColor = colorCache.get(clusterGroup);
      }

      let foregroundColor = "rgb(0, 0, 0)";
      if (!colorCache.has(`${backgroundColor}-fg`)) {
        foregroundColor = tinycolor(backgroundColor).isLight()
          ? "#000"
          : "#fff";
        colorCache.set(`${backgroundColor}-fg`, foregroundColor);
      } else {
        foregroundColor = colorCache.get(`${backgroundColor}-fg`);
      }

      const roundedValue = (value as number).toFixed(2);

      return {
        x,
        y,
        value: value as number,
        displayValue: value === 100 ? "100" : roundedValue.toString(),
        backgroundColor,
        foregroundColor,
      };
    }),
  );

  return result;
};
