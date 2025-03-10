import tinycolor from "tinycolor2";
import type { ColorScaleArray } from "./colorScales";
import { createD3ColorScale, distinctColor } from "./colors";
import type { HeatmapData, HeatmapSettings } from "./plotTypes";

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
  clusterData?: { id: string; group: number }[],
) =>
  data.flatMap((row, y) =>
    row.filter(Number).map((value, x) => {
      const clusterX = clusterData?.find((i) => i.id === tickText[x])?.group;
      const clusterY = clusterData?.find((i) => i.id === tickText[y])?.group;

      const clusterMatch =
        clusterX !== undefined &&
        clusterY !== undefined &&
        clusterX === clusterY;

      const clusterGroup = clusterMatch ? clusterX : null;

      const backgroundColor = clusterGroup
        ? distinctColor(clusterGroup)
        : "rgb(245, 245, 245)";
      const foregroundColor = tinycolor(backgroundColor).isLight()
        ? "#000"
        : "#fff";
      const roundedValue = Number(Number(value).toFixed(2));

      return {
        x,
        y,
        value: roundedValue,
        displayValue: roundedValue.toString(),
        backgroundColor,
        foregroundColor,
      };
    }),
  );
