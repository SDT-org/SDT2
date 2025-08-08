import tinycolor from "tinycolor2";
import type { ColorScaleArray } from "../../../colorScales";
import { createD3ColorScale, distinctColor } from "../../../colors";
import type {
  ClusterDataItem,
  HeatmapData,
  HeatmapSettings,
} from "../../../plotTypes";

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
    | "colorScaleKey"
    | "vmax"
    | "vmin"
    | "annotation_rounding"
    | "hideValuesBelow"
    | "hideValuesBelowEnabled"
  >,
  colorScale: ColorScaleArray,
) => {
  const colorFn = createD3ColorScale(
    colorScale,
    settings.colorScaleKey === "Discrete",
    settings.vmax,
    settings.vmin,
  );

  const emptyColor = "#f5f5f5";
  const blackColor = "#000";
  const whiteColor = "#fff";
  const hideValuesBelowEnabled = settings.hideValuesBelowEnabled;
  const hideValuesBelow = settings.hideValuesBelow;

  const colorCache = new Map<string, { bg: string; fg: string }>();

  // Pre-allocate result array to avoid multiple array resizing
  let totalCells = 0;
  for (let y = 0; y < data.length; y++) {
    const row = data[y];
    if (!row) continue;
    for (let x = 0; x < row.length; x++) {
      if (row[x] === 0 || Number(row[x])) {
        totalCells++;
      }
    }
  }

  const result = new Array(totalCells);
  let resultIndex = 0;

  for (let y = 0; y < data.length; y++) {
    const row = data[y];
    if (!row) continue;
    let actualX = 0;

    for (let x = 0; x < row.length; x++) {
      const datum = row[x];

      // can skip empty/null values early
      if (!(datum === 0 || Number(datum))) {
        continue;
      }

      const value = datum as number;

      const shouldHide = hideValuesBelowEnabled && value < hideValuesBelow;
      const effectiveValue = shouldHide ? 0 : value;

      let displayValue: string;
      if (effectiveValue === 100) {
        displayValue = "100";
      } else if (effectiveValue === 0) {
        displayValue = "";
      } else {
        displayValue = effectiveValue.toFixed(settings.annotation_rounding);
      }

      let backgroundColor: string;
      let foregroundColor: string;

      if (displayValue === "") {
        backgroundColor = emptyColor;
        foregroundColor = blackColor;
      } else {
        const colorKey = effectiveValue.toString();
        const colors = colorCache.get(colorKey);

        if (!colors) {
          backgroundColor = colorFn(effectiveValue);
          foregroundColor = tinycolor(backgroundColor).isLight()
            ? blackColor
            : whiteColor;
          colorCache.set(colorKey, {
            bg: backgroundColor,
            fg: foregroundColor,
          });
        } else {
          backgroundColor = colors.bg;
          foregroundColor = colors.fg;
        }
      }

      result[resultIndex] = {
        x: actualX,
        y,
        value,
        displayValue,
        backgroundColor,
        foregroundColor,
      };

      resultIndex++;
      actualX++;
    }
  }

  return result;
};

export const formatClustermapData = (
  data: HeatmapData,
  tickText: string[],
  clusterData?: ClusterDataItem[],
) => {
  const clusterMap = new Map();
  const originalClusterMap = new Map();
  if (clusterData) {
    for (const item of clusterData) {
      clusterMap.set(item.id, item.cluster);
      // Use original_cluster if available, otherwise use cluster
      originalClusterMap.set(item.id, item.original_cluster ?? item.cluster);
    }
  }

  const colorCache = new Map();

  const result = data.flatMap((row, y) =>
    row
      .filter((datum) => datum === 0 || Number(datum))
      .map((value, x) => {
        const clusterX = clusterMap.get(tickText[x]);
        const clusterY = clusterMap.get(tickText[y]);
        const originalClusterX = originalClusterMap.get(tickText[x]);
        const originalClusterY = originalClusterMap.get(tickText[y]);

        const clusterMatch =
          clusterX !== undefined &&
          clusterY !== undefined &&
          clusterX === clusterY;

        const clusterGroup = clusterMatch ? clusterX : null;
        const originalClusterGroup =
          clusterMatch && originalClusterX === originalClusterY
            ? originalClusterX
            : null;

        let backgroundColor = "rgb(245, 245, 245)";
        if (clusterGroup && originalClusterGroup) {
          if (!colorCache.has(clusterGroup)) {
            // Use original cluster number for color generation
            colorCache.set(clusterGroup, distinctColor(originalClusterGroup));
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
        const displayValue =
          value === 100 ? "100" : value === 0 ? "" : roundedValue.toString();

        return {
          x,
          y,
          value: value as number,
          displayValue,
          backgroundColor,
          foregroundColor,
        };
      }),
  );

  return result;
};
