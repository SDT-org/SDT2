import * as d3 from "d3";
import { distinctColor } from "../../../colors";
import { plotFontMonospace } from "../../../constants";
import type { HeatmapRenderProps } from "./Heatmap";
import { getCellMetrics } from "./heatmapUtils";

interface CanvasRenderParams
  extends Omit<
    HeatmapRenderProps,
    | "onRenderComplete"
    | "showLegend"
    | "showClusterCounts"
    | "clusterData"
    | "clusterCounts"
  > {
  canvas: HTMLCanvasElement;
  transform: d3.ZoomTransform;
  showLegend: HeatmapRenderProps["showLegend"] | undefined;
  showClusterCounts: HeatmapRenderProps["showClusterCounts"] | undefined;
  clusterData: HeatmapRenderProps["clusterData"] | undefined;
  clusterCounts: HeatmapRenderProps["clusterCounts"] | undefined;
}

export const renderHeatmapCanvas = ({
  canvas,
  data,
  clusterData,
  tickText,
  colorScale,
  minVal,
  maxVal,
  width,
  height,
  cellSpace,
  roundTo,
  cbarHeight,
  cbarWidth,
  axlabel_fontsize,
  axlabel_xrotation,
  axlabel_yrotation,
  titleFont,
  showPercentIdentities,
  showTitles,
  title,
  axis_labels,
  showscale,
  margin,
  settings,
  showLegend,
  showClusterCounts,
  clusterCounts,
  transform,
}: CanvasRenderParams): void => {
  const ctx = canvas.getContext("2d");
  if (!ctx) {
    console.warn("Failed to get 2d context");
    return;
  }

  const pixelRatio = window.devicePixelRatio || 1;
  canvas.width = width * pixelRatio;
  canvas.height = height * pixelRatio;
  canvas.style.width = `${width}px`;
  canvas.style.height = `${height}px`;

  ctx.fillStyle = "#fff";
  ctx.fillRect(0, 0, canvas.width, canvas.height);
  ctx.scale(pixelRatio, pixelRatio);

  ctx.textRendering = "optimizeSpeed";

  const plotSize = Math.min(width, height) - margin.left - margin.right;
  const cellSize = plotSize / tickText.length;
  const scale = d3
    .scaleLinear()
    .domain([maxVal, minVal])
    .range([0, cbarHeight]);
  const tickValues = scale.ticks(5);

  // set zoom transform
  ctx.save();
  ctx.translate(transform.x + margin.left, transform.y + margin.top);
  ctx.scale(transform.k, transform.k);

  const cellMetrics = getCellMetrics(cellSize, cellSpace, roundTo + 3);
  ctx.textAlign = "center";
  ctx.textBaseline = "middle";
  ctx.font = `${cellMetrics.fontSize}px ${plotFontMonospace.family}`;

  //index data
  const rows = new Map();
  const cols = new Map();
  let rowIndex = 0;
  let colIndex = 0;

  for (const { x, y } of data) {
    if (!rows.has(x)) rows.set(x, rowIndex++);
    if (!cols.has(y)) cols.set(y, colIndex++);
  }

  // Draw cells
  for (const d of data) {
    const x = rows.get(d.x) * cellSize + cellMetrics.cellOffset;
    const y = cols.get(d.y) * cellSize + cellMetrics.cellOffset;

    ctx.fillStyle = d.backgroundColor;
    ctx.fillRect(x, y, cellMetrics.cellSize, cellMetrics.cellSize);

    if (showPercentIdentities) {
      ctx.fillStyle = d.foregroundColor;
      ctx.fillText(
        d.displayValue,
        x + cellMetrics.textOffset,
        y + cellMetrics.textOffset,
      );
    }
  }
  ctx.restore();

  if (showTitles) {
    ctx.fillStyle = "black";
    ctx.textAlign = "center";
    ctx.textBaseline = "top";
    ctx.font = `Bold 20px ${titleFont.family}`;
    ctx.fillText(title, width / 2, margin.top - 20);
  }

  if (axis_labels) {
    const axisGap = 5;

    // X-axis labels
    for (const [i, txt] of tickText.entries()) {
      if (txt === undefined) continue;
      ctx.save();
      ctx.translate(
        margin.left + // X starting position
          i * cellSize * transform.k + // X position for current tick (accounting for zoom)
          (cellSize * transform.k) / 2 + // Center tick horizontally within cell
          transform.x, // X pan offset for cell
        transform.y + margin.top + plotSize * transform.k + axisGap,
      );
      ctx.rotate(((axlabel_xrotation + 270) * Math.PI) / 180);
      ctx.fillStyle = "black";
      ctx.textAlign = "right";
      ctx.textBaseline = "middle";
      ctx.font = `${Math.max(
        axlabel_fontsize * transform.k,
        axlabel_fontsize,
      )}px ${plotFontMonospace.family}`;
      ctx.fillText(txt, 0, 0);
      ctx.restore();
    }

    // Y-axis labels
    for (const [i, txt] of tickText.entries()) {
      if (txt === undefined) continue;
      ctx.save();
      ctx.translate(
        transform.x + margin.left - axisGap,
        margin.top + // Y start position (top margin)
          i * cellSize * transform.k + // Y position for current tick (accounting for zoom)
          (cellSize * transform.k) / 2 + // Center vertically within cell
          transform.y, // Y pan offset
      );
      ctx.rotate(((axlabel_yrotation + 360) * Math.PI) / 180);
      ctx.fillStyle = "black";
      ctx.textAlign = "right";
      ctx.textBaseline = "middle";
      ctx.font = `${Math.max(
        axlabel_fontsize * transform.k,
        axlabel_fontsize,
      )}px ${plotFontMonospace.family}`;
      ctx.fillText(txt, 0, 0);
      ctx.restore();
    }
  }

  if (showscale) {
    // sticky the colorbar to the right side of the canvas or plot, depending on the space available
    const positionX = Math.min(
      margin.left + plotSize + cbarWidth,
      width - cbarWidth - margin.right,
    );
    const gradient = ctx.createLinearGradient(
      margin.left,
      margin.top + cbarHeight,
      margin.left,
      margin.top,
    );
    // addded normalize to min max to fix gradient not matching colorscale
    for (const [stop, color] of colorScale) {
      let stopValue = stop;
      if (settings?.colorScaleKey === "Discrete") {
        stopValue = (stop - minVal) / (maxVal - minVal);
      }
      gradient.addColorStop(stopValue, color);
    }

    ctx.fillStyle = gradient;
    ctx.fillRect(positionX, margin.top, cbarWidth, cbarHeight);

    ctx.fillStyle = "black";
    ctx.textAlign = "left";
    ctx.textBaseline = "middle";
    ctx.font = "10px 'Roboto Mono'";

    for (const tick of tickValues) {
      const y = scale(tick) + margin.top;
      ctx.fillText(tick.toString(), positionX + cbarWidth + 9, y);
      ctx.beginPath();
      ctx.moveTo(positionX + cbarWidth, y);
      ctx.lineTo(positionX + cbarWidth + 6, y);
      ctx.stroke();
    }
  }

  if (clusterData && showLegend) {
    const legendWidth = showClusterCounts ? 120 : 80;
    const cellSize = 10;
    const lineGap = 20;
    const labelGap = 5;
    const columnGap = showClusterCounts ? 30 : 20;

    // sticky the legend to the right side of the canvas or plot, depending on the space available
    const positionX = Math.min(
      margin.left + plotSize,
      width - legendWidth * 2 - columnGap - margin.right,
    );

    const uniqueClusters = [...new Set(clusterData.map((i) => i.cluster))]
      .sort((a, b) => a - b)
      .slice(0, 50); // Arbitrarily set to 50

    // Create a map of cluster to original cluster for color consistency
    const clusterToOriginal = new Map();
    for (const item of clusterData) {
      if (!clusterToOriginal.has(item.cluster)) {
        clusterToOriginal.set(
          item.cluster,
          item.original_cluster ?? item.cluster,
        );
      }
    }

    const clusterColors = uniqueClusters.map((cluster) =>
      distinctColor(clusterToOriginal.get(cluster) ?? cluster),
    );
    ctx.font = `10px 'Roboto Mono'`;
    ctx.textAlign = "left";
    ctx.textBaseline = "middle";

    uniqueClusters.forEach((cluster, index) => {
      // Determine column (0  left, 1  right)
      const column = index % 2;

      const row = Math.floor(index / 2);
      const textX = positionX + column * (legendWidth + columnGap);
      const textY = margin.top + lineGap * row;

      ctx.fillStyle = clusterColors[index] as string;
      ctx.fillRect(textX, textY, cellSize, cellSize);

      ctx.fillStyle = "black";
      // Format cluster number with leading space for single digits
      const clusterStr = cluster < 10 ? ` ${cluster}` : `${cluster}`;

      if (showClusterCounts && clusterCounts?.[cluster]) {
        // Add leading space before bracket for single digits
        const countDisplay =
          clusterCounts[cluster] < 10
            ? ` [${clusterCounts[cluster]}]`
            : `[${clusterCounts[cluster]}]`;

        // Draw cluster label
        ctx.textAlign = "left";
        ctx.fillText(
          `Cluster ${clusterStr}`,
          textX + cellSize + labelGap,
          textY + cellSize / 2,
        );
        // Draw count right-aligned
        ctx.textAlign = "right";
        ctx.fillText(countDisplay, textX + legendWidth, textY + cellSize / 2);
        ctx.textAlign = "left"; // Reset alignment
      } else {
        ctx.fillText(
          `Cluster ${clusterStr}`,
          textX + cellSize + labelGap,
          textY + cellSize / 2,
        );
      }
    });
  }
};
