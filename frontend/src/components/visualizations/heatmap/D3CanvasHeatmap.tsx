import * as d3 from "d3";
import React from "react";
import { useExportCanvas } from "../../../hooks/useExportCanvas";
import type { MetaData } from "../../../plotTypes";
import type { HeatmapRenderProps } from "./Heatmap";
import { renderHeatmapCanvas } from "./canvasRenderer";

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

  const drawing = React.useRef(false);

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

    renderHeatmapCanvas({
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
    });

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
    margin,
    minVal,
    maxVal,
    settings,
    clusterData,
    showLegend,
    showClusterCounts,
    clusterCounts,
    onRenderComplete,
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
  }, [width, height, margin]);

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
      clusterData.find((i) => i.id === tickText[cell.x])?.cluster ===
        clusterData.find((i) => i.id === tickText[cell.y])?.cluster
        ? clusterData.find((i) => i.id === tickText[cell.x])?.cluster
        : null;

    if (cell) {
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
  };

  const idValue =
    tooltipData &&
    (clusterData ? tooltipData.percentId : tooltipData.value)?.toFixed(2);

  const idDisplay = Number(idValue) === 0 ? "" : idValue;

  return (
    <div style={{ position: "relative" }}>
      <canvas
        id="heatmap-canvas"
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
