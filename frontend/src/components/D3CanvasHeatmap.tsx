import * as d3 from "d3";
import React from "react";
import { distinctColor } from "../colors";
import { plotFontMonospace } from "../constants";
import { getCellMetrics } from "../heatmapUtils";
import { useHeatmapRef } from "../hooks/useHeatmapRef";
import type { HeatmapRenderProps } from "./Heatmap";

export const D3CanvasHeatmap = ({
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
}: HeatmapRenderProps) => {
  const canvasRef =
    useHeatmapRef() as React.MutableRefObject<HTMLCanvasElement>;
  const [transform, setTransform] = React.useState(d3.zoomIdentity);
  const [tooltipData, setTooltipData] = React.useState<{
    x: number;
    y: number;
    value: number | null;
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

  const drawCanvas = React.useCallback(() => {
    const canvas = canvasRef.current;
    if (!canvas) {
      console.warn("Failed to find canvas");
      return;
    }

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

    ctx.clearRect(0, 0, canvas.width, canvas.height);
    ctx.scale(pixelRatio, pixelRatio);

    ctx.textRendering = "optimizeSpeed";

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
      const positionX = width - cbarWidth - margin.right;
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
      const legendWidth = 80;
      const cellSize = 10;
      const lineGap = 20;
      const labelGap = 5;
      const columnGap = 20;
      const positionX = width - legendWidth * 2 - columnGap - margin.right;

      const uniqueClusters = [
        ...new Set(clusterData.map((i) => i.group)),
      ].slice(0, 50);

      uniqueClusters.forEach((cluster, index) => {
        // Determine column (0 for left, 1 for right)
        const column = index % 2;

        // Calculate row position (every two items share the same row)
        const row = Math.floor(index / 2);

        // Calculate position based on column and row
        const itemX = positionX + column * (legendWidth + columnGap);
        const itemY = margin.top + lineGap * row;

        // Draw colored square
        ctx.fillStyle = distinctColor(index + 1);
        ctx.fillRect(itemX, itemY, cellSize, cellSize);

        // Draw text
        ctx.fillStyle = "black";
        ctx.textAlign = "left";
        ctx.textBaseline = "middle";
        ctx.font = `10px 'Roboto Mono'`;
        ctx.fillText(
          `Cluster ${cluster.toString()}`,
          itemX + cellSize + labelGap,
          itemY + cellSize / 2,
        );
      });
    }
  }, [
    transform,
    data,
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
    canvasRef,
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
  ]);

  React.useEffect(() => {
    drawCanvas();
  }, [drawCanvas]);

  React.useEffect(() => {
    const canvas = canvasRef.current;
    if (!canvas) return;

    const zoom = d3
      .zoom()
      .scaleExtent([0.5, 10])
      .translateExtent([
        [-margin.left, -margin.top],
        [width, height],
      ])
      .on("zoom", (event) => setTransform(event.transform));

    d3.select(canvas).call(
      zoom as unknown as d3.ZoomBehavior<HTMLCanvasElement, unknown>,
    );

    return () => {
      d3.select(canvas).on("zoom", null);
    };
  }, [canvasRef, width, height, margin]);

  const handleMouseMove = (event: React.MouseEvent<HTMLCanvasElement>) => {
    const canvas = canvasRef.current;
    if (!canvas) return;

    const rect = canvas.getBoundingClientRect();
    const x = event.clientX - rect.left;
    const y = event.clientY - rect.top;

    const dataX = Math.floor(
      (x - margin.left - transform.x) / (cellSize * transform.k),
    );
    const dataY = Math.floor(
      (y - margin.top - transform.y) / (cellSize * transform.k),
    );

    const cell = data.find((d) => d.x === dataX && d.y === dataY);

    const clusterGroup =
      clusterData &&
      cell &&
      clusterData.find((i) => i.id === tickText[cell.x])?.group ===
        clusterData.find((i) => i.id === tickText[cell.y])?.group
        ? clusterData.find((i) => i.id === tickText[cell.x])?.group
        : null;

    if (cell) {
      setTooltipData({
        x,
        y,
        value: clusterData ? (clusterGroup ?? null) : cell.value,
        xLabel: tickText[cell.x] || "",
        yLabel: tickText[cell.y] || "",
      });
    } else {
      setTooltipData(null);
    }
  };

  return (
    <div style={{ position: "relative" }}>
      <canvas
        ref={canvasRef}
        style={{ background: "#fff" }}
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
            <dt>SeqX:</dt>
            <dd>{tooltipData.xLabel}</dd>
          </div>
          <div>
            <dt> SeqY:</dt>
            <dd>{tooltipData.yLabel}</dd>
          </div>
          <div>
            <dt>
              {clusterData
                ? tooltipData.value
                  ? "Cluster:"
                  : ""
                : "Percent ID:"}
            </dt>
            <dd>
              {clusterData
                ? tooltipData.value || ""
                : `${tooltipData?.value?.toFixed(2)}%`}
            </dd>
          </div>
        </dl>
      )}
    </div>
  );
};
