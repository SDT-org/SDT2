import * as d3 from "d3";
import React from "react";
import { distinctColor } from "../colors";
import { plotFontMonospace } from "../constants";
import { getCellMetrics } from "../heatmapUtils";
import { useExportCanvas } from "../hooks/useExportCanvas";
import type { MetaData } from "../plotTypes";
import type { HeatmapRenderProps } from "./Heatmap";

const getMetricLabel = (metaData?: MetaData): string => {
  if (!metaData?.run) {
    return "Identity";
  }

  if (
    metaData.run.analysis_method === "lzani" &&
    metaData.run.lzani?.score_type
  ) {
    const scoreType = metaData.run.lzani.score_type.toUpperCase();
    return scoreType === "ANI"
      ? "ANI"
      : scoreType === "TANI"
        ? "TANI"
        : "Identity";
  }

  return "Percent ID";
};

// Simple debounce function
function debounce<T extends (...args: unknown[]) => unknown>(
  func: T,
  wait: number,
): (...args: Parameters<T>) => void {
  let timeout: number | null = null;
  return (...args: Parameters<T>) => {
    if (timeout) clearTimeout(timeout);
    timeout = window.setTimeout(() => func(...args), wait);
  };
}

export const D3CanvasHeatmapOptimized = ({
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
  onRenderComplete,
  metaData,
}: HeatmapRenderProps) => {
  const canvasRef = React.useRef<HTMLCanvasElement>(null);
  const exportCanvas = useExportCanvas(
    clusterData ? "clustermap" : "heatmap",
    canvasRef,
  );

  const [transform, setTransform] = React.useState(d3.zoomIdentity);
  const [tooltipData, setTooltipData] = React.useState<{
    x: number;
    y: number;
    value: number | null;
    percentId?: number | null;
    cluster?: number | null;
    xLabel?: string;
    yLabel?: string;
  } | null>(null);

  const plotSize = Math.min(width, height) - margin.left - margin.right;
  const cellSize = plotSize / tickText.length;

  const scale = React.useMemo(
    () => d3.scaleLinear().domain([maxVal, minVal]).range([0, cbarHeight]),
    [minVal, maxVal, cbarHeight],
  );

  const tickValues = React.useMemo(() => scale.ticks(5), [scale]);

  const drawing = React.useRef(false);
  const rafId = React.useRef<number | null>(null);

  // Pre-index data for faster lookups
  const dataIndex = React.useMemo(() => {
    const index = new Map<string, (typeof data)[0]>();
    for (const d of data) {
      index.set(`${d.x},${d.y}`, d);
    }
    return index;
  }, [data]);

  // Calculate visible bounds for chunked rendering
  const calculateVisibleBounds = React.useCallback(
    (transform: d3.ZoomTransform) => {
      const bufferCells = 50; // Buffer around viewport
      const effectiveCellSize = cellSize * transform.k;

      return {
        startX: Math.max(
          0,
          Math.floor((-transform.x - margin.left) / effectiveCellSize) -
            bufferCells,
        ),
        endX: Math.min(
          tickText.length,
          Math.ceil((-transform.x + width - margin.left) / effectiveCellSize) +
            bufferCells,
        ),
        startY: Math.max(
          0,
          Math.floor((-transform.y - margin.top) / effectiveCellSize) -
            bufferCells,
        ),
        endY: Math.min(
          tickText.length,
          Math.ceil((-transform.y + height - margin.top) / effectiveCellSize) +
            bufferCells,
        ),
      };
    },
    [cellSize, width, height, margin, tickText.length],
  );

  const drawCanvas = React.useCallback(() => {
    if (drawing.current) {
      console.warn("Drawing already in progress");
      return;
    }
    drawing.current = true;

    const canvas = canvasRef.current;
    if (!canvas) {
      console.warn("Failed to find canvas");
      drawing.current = false;
      return;
    }

    const ctx = canvas.getContext("2d", { alpha: false });
    if (!ctx) {
      console.warn("Failed to get 2d context");
      drawing.current = false;
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

    // Disable anti-aliasing for sharper rendering
    ctx.imageSmoothingEnabled = false;

    // set zoom transform
    ctx.save();
    ctx.translate(transform.x + margin.left, transform.y + margin.top);
    ctx.scale(transform.k, transform.k);

    const cellMetrics = getCellMetrics(cellSize, cellSpace, roundTo + 3);

    // Calculate visible bounds for chunked rendering
    const bounds = calculateVisibleBounds(transform);

    // Only render visible cells
    const visibleCellCount =
      (bounds.endX - bounds.startX) * (bounds.endY - bounds.startY);
    const useChunkedRendering = visibleCellCount < data.length * 0.5; // Use chunking if rendering less than 50% of cells

    if (useChunkedRendering) {
      // Chunked rendering for zoomed-in views
      ctx.textAlign = "center";
      ctx.textBaseline = "middle";
      ctx.font = `${cellMetrics.fontSize}px ${plotFontMonospace.family}`;

      for (let x = bounds.startX; x < bounds.endX; x++) {
        for (let y = bounds.startY; y < bounds.endY; y++) {
          const d = dataIndex.get(`${x},${y}`);
          if (!d) continue;

          const xPos = x * cellSize + cellMetrics.cellOffset;
          const yPos = y * cellSize + cellMetrics.cellOffset;

          ctx.fillStyle = d.backgroundColor;
          ctx.fillRect(xPos, yPos, cellMetrics.cellSize, cellMetrics.cellSize);

          if (showPercentIdentities && cellMetrics.fontSize > 2) {
            ctx.fillStyle = d.foregroundColor;
            ctx.fillText(
              d.displayValue,
              xPos + cellMetrics.textOffset,
              yPos + cellMetrics.textOffset,
            );
          }
        }
      }
    } else {
      // Batch rendering for zoomed-out views
      // Group cells by color for batch rendering
      const colorGroups = new Map<string, Array<{ x: number; y: number }>>();

      for (const d of data) {
        if (!colorGroups.has(d.backgroundColor)) {
          colorGroups.set(d.backgroundColor, []);
        }
        const group = colorGroups.get(d.backgroundColor);
        if (group) {
          group.push({
            x: d.x * cellSize + cellMetrics.cellOffset,
            y: d.y * cellSize + cellMetrics.cellOffset,
          });
        }
      }

      // Render all cells of the same color at once
      for (const [color, positions] of colorGroups) {
        ctx.fillStyle = color;
        for (const pos of positions) {
          ctx.fillRect(
            pos.x,
            pos.y,
            cellMetrics.cellSize,
            cellMetrics.cellSize,
          );
        }
      }

      // Only render text if cells are large enough
      if (showPercentIdentities && cellMetrics.fontSize > 2) {
        ctx.textAlign = "center";
        ctx.textBaseline = "middle";
        ctx.font = `${cellMetrics.fontSize}px ${plotFontMonospace.family}`;

        for (const d of data) {
          if (d.displayValue) {
            ctx.fillStyle = d.foregroundColor;
            ctx.fillText(
              d.displayValue,
              d.x * cellSize + cellMetrics.cellOffset + cellMetrics.textOffset,
              d.y * cellSize + cellMetrics.cellOffset + cellMetrics.textOffset,
            );
          }
        }
      }
    }

    ctx.restore();

    // Draw static elements (title, axes, scale, legend)
    if (showTitles) {
      ctx.fillStyle = "black";
      ctx.textAlign = "center";
      ctx.textBaseline = "top";
      ctx.font = `Bold 20px ${titleFont.family}`;
      ctx.fillText(title, width / 2, margin.top - 20);
    }

    if (axis_labels) {
      const axisGap = 5;

      // Only render visible axis labels
      const labelBounds = calculateVisibleBounds(transform);

      // X-axis labels
      for (let i = labelBounds.startX; i < labelBounds.endX; i++) {
        const txt = tickText[i];
        if (txt === undefined) continue;
        ctx.save();
        ctx.translate(
          margin.left +
            i * cellSize * transform.k +
            (cellSize * transform.k) / 2 +
            transform.x,
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
      for (let i = labelBounds.startY; i < labelBounds.endY; i++) {
        const txt = tickText[i];
        if (txt === undefined) continue;
        ctx.save();
        ctx.translate(
          transform.x + margin.left - axisGap,
          margin.top +
            i * cellSize * transform.k +
            (cellSize * transform.k) / 2 +
            transform.y,
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

      const positionX = Math.min(
        margin.left + plotSize,
        width - legendWidth * 2 - columnGap - margin.right,
      );

      const uniqueClusters = [...new Set(clusterData.map((i) => i.cluster))]
        .sort((a, b) => a - b)
        .slice(0, 50);

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
        const column = index % 2;
        const row = Math.floor(index / 2);
        const textX = positionX + column * (legendWidth + columnGap);
        const textY = margin.top + lineGap * row;

        ctx.fillStyle = clusterColors[index] as string;
        ctx.fillRect(textX, textY, cellSize, cellSize);

        ctx.fillStyle = "black";
        const clusterStr = cluster < 10 ? ` ${cluster}` : `${cluster}`;

        if (showClusterCounts && clusterCounts?.[cluster]) {
          const countDisplay =
            clusterCounts[cluster] < 10
              ? ` [${clusterCounts[cluster]}]`
              : `[${clusterCounts[cluster]}]`;

          ctx.textAlign = "left";
          ctx.fillText(
            `Cluster ${clusterStr}`,
            textX + cellSize + labelGap,
            textY + cellSize / 2,
          );
          ctx.textAlign = "right";
          ctx.fillText(countDisplay, textX + legendWidth, textY + cellSize / 2);
          ctx.textAlign = "left";
        } else {
          ctx.fillText(
            `Cluster ${clusterStr}`,
            textX + cellSize + labelGap,
            textY + cellSize / 2,
          );
        }
      });
    }

    drawing.current = false;

    const id = setTimeout(() => {
      exportCanvas();
    }, 0);

    onRenderComplete?.();

    return () => {
      clearTimeout(id);
    };
  }, [
    exportCanvas,
    transform,
    data,
    dataIndex,
    scale,
    tickValues,
    cellSize,
    cellSpace,
    showPercentIdentities,
    roundTo,
    showTitles,
    title,
    axlabel_fontsize,
    axlabel_xrotation,
    axlabel_yrotation,
    titleFont,
    tickText,
    cbarWidth,
    cbarHeight,
    colorScale,
    width,
    height,
    axis_labels,
    showscale,
    plotSize,
    margin,
    minVal,
    maxVal,
    settings?.colorScaleKey,
    clusterData,
    showLegend,
    showClusterCounts,
    clusterCounts,
    onRenderComplete,
    calculateVisibleBounds,
  ]);

  // Debounced render function
  const debouncedRender = React.useMemo(
    () =>
      debounce(() => {
        if (rafId.current) {
          cancelAnimationFrame(rafId.current);
        }
        rafId.current = requestAnimationFrame(() => {
          drawCanvas();
          rafId.current = null;
        });
      }, 16), // ~60fps
    [drawCanvas],
  );

  React.useEffect(() => {
    debouncedRender();
    return () => {
      if (rafId.current) {
        cancelAnimationFrame(rafId.current);
      }
    };
  }, [debouncedRender]);

  React.useEffect(() => {
    const canvas = canvasRef.current;
    if (!canvas) return;

    const zoom = d3
      .zoom()
      .scaleExtent([0.1, 20]) // Allow more zoom out for large datasets
      .translateExtent([
        [-margin.left * 2, -margin.top * 2],
        [width + margin.right * 2, height + margin.bottom * 2],
      ])
      .on("zoom", (event) => {
        setTransform(event.transform);
      });

    d3.select(canvas).call(
      zoom as unknown as d3.ZoomBehavior<HTMLCanvasElement, unknown>,
    );

    return () => {
      d3.select(canvas).on("zoom", null);
    };
  }, [width, height, margin]);

  // Optimized mouse move handler with RAF throttling
  const handleMouseMove = React.useMemo(() => {
    let rafId: number | null = null;
    let lastEvent: React.MouseEvent<HTMLCanvasElement> | null = null;

    const processMouseMove = () => {
      if (!lastEvent) return;

      const canvas = canvasRef.current;
      if (!canvas) return;

      const rect = canvas.getBoundingClientRect();
      const x = lastEvent.clientX - rect.left;
      const y = lastEvent.clientY - rect.top;

      const dataX = Math.floor(
        (x - margin.left - transform.x) / (cellSize * transform.k),
      );
      const dataY = Math.floor(
        (y - margin.top - transform.y) / (cellSize * transform.k),
      );

      // Use indexed lookup instead of array search
      const cell = dataIndex.get(`${dataX},${dataY}`);

      if (cell) {
        const clusterGroup =
          clusterData &&
          clusterData.find((i) => i.id === tickText[cell.x])?.cluster ===
            clusterData.find((i) => i.id === tickText[cell.y])?.cluster
            ? clusterData.find((i) => i.id === tickText[cell.x])?.cluster
            : null;

        setTooltipData({
          x,
          y,
          value: clusterData ? (clusterGroup ?? null) : cell.value,
          percentId: cell.value,
          cluster: clusterGroup ?? null,
          xLabel: tickText[cell.x] || "",
          yLabel: tickText[cell.y] || "",
        });
      } else {
        setTooltipData(null);
      }

      rafId = null;
    };

    return (event: React.MouseEvent<HTMLCanvasElement>) => {
      lastEvent = event;
      if (!rafId) {
        rafId = requestAnimationFrame(processMouseMove);
      }
    };
  }, [cellSize, transform, margin, dataIndex, tickText, clusterData]);

  const idValue =
    tooltipData &&
    (clusterData ? tooltipData.percentId : tooltipData.value)?.toFixed(2);

  const idDisplay = Number(idValue) === 0 ? "" : idValue;

  return (
    <div style={{ position: "relative" }}>
      <canvas
        id="heatmap-canvas"
        ref={canvasRef}
        style={{ background: "#fff", cursor: "crosshair" }}
        width={width}
        height={height}
        onMouseMove={handleMouseMove}
        onMouseLeave={() => setTooltipData(null)}
      />
      {tooltipData && (
        <dl
          className="heatmap-tooltip"
          style={{
            left: tooltipData.x + 10,
            top: tooltipData.y + 10,
          }}
        >
          <div>
            <dt>Seq X:</dt>
            <dd>{tooltipData.xLabel}</dd>
          </div>
          <div>
            <dt>Seq Y:</dt>
            <dd>{tooltipData.yLabel}</dd>
          </div>
          {clusterData ? (
            <>
              {tooltipData.cluster !== null && (
                <div>
                  <dt>Cluster:</dt>
                  <dd>{tooltipData.cluster}</dd>
                </div>
              )}
            </>
          ) : null}
          <div>
            <dt>{getMetricLabel(metaData)}:</dt>
            <dd>{idDisplay ? `${idDisplay}%` : "Unaligned"}</dd>
          </div>
        </dl>
      )}
    </div>
  );
};
