import * as d3 from "d3";
import React from "react";
import tinycolor from "tinycolor2";
import type { ColorScaleArray } from "../colorScales";
import { plotFontMonospace } from "../constants";
import { useHeatmapRef } from "../hooks/useHeatmapRef";
import type { HeatmapRenderProps } from "./Heatmap";

function createD3ColorScale(
  colorArray: ColorScaleArray,
  minValue: number,
  maxValue: number,
): d3.ScaleLinear<string, string> {
  const domain = colorArray.map(
    ([stop]) => stop * (maxValue - minValue) + minValue,
  );
  const range = colorArray.map(([_, color]) => color);

  return d3
    .scaleLinear<string>()
    .domain(domain)
    .range(range)
    .interpolate(d3.interpolateRgb);
}

export const D3CanvasHeatmap = ({
  data,
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
  annotation_font_size,
  axlabel_xfontsize,
  axlabel_yfontsize,
  axlabel_xrotation,
  axlabel_yrotation,
  titleFont,
  showPercentIdentities,
  showTitles,
  title,
  subtitle,
  axis_labels,
  showscale,
  margin,
}: HeatmapRenderProps) => {
  const canvasRef =
    useHeatmapRef() as React.MutableRefObject<HTMLCanvasElement>;
  const [transform, setTransform] = React.useState(d3.zoomIdentity);
  const [tooltipData, setTooltipData] = React.useState<{
    x: number;
    y: number;
    value: number;
    xLabel?: string;
    yLabel?: string;
  } | null>(null);

  const filteredData = React.useMemo(
    () => data.filter((d) => Number(d.value)),
    [data],
  );

  const size = Math.min(width, height);
  const plotSize = size - margin.left - margin.right;

  const n = tickText.length;
  const cellSize = plotSize / n;

  const colorFn = React.useMemo(
    () => createD3ColorScale(colorScale, minVal, maxVal),
    [colorScale, minVal, maxVal],
  );

  const scale = React.useMemo(
    () => d3.scaleLinear().domain([maxVal, minVal]).range([0, cbarHeight]),
    [minVal, maxVal, cbarHeight],
  );

  const tickValues = React.useMemo(() => scale.ticks(5), [scale]);
  const minValue = Math.min(...filteredData.map((item) => item.value));

  const drawCanvas = React.useCallback(() => {
    // console.count("drawCanvas");

    const canvas = canvasRef.current;
    if (!canvas) return;

    const ctx = canvas.getContext("2d");
    if (!ctx) return;

    const pixelRatio = window.devicePixelRatio || 1;
    canvas.width = width * pixelRatio;
    canvas.height = height * pixelRatio;
    canvas.style.width = `${width}px`;
    canvas.style.height = `${height}px`;

    ctx.clearRect(0, 0, canvas.width, canvas.height);
    ctx.scale(pixelRatio, pixelRatio);

    // set zoom transform
    ctx.save();
    ctx.translate(transform.x + margin.left, transform.y + margin.top);
    ctx.scale(transform.k, transform.k);
    //indexz data
    const rows = [...new Set(filteredData.map((d) => d.x))];
    const cols = [...new Set(filteredData.map((d) => d.y))];

    // Draw cells
    for (const d of filteredData) {
      const x = cols.indexOf(d.x) * cellSize + cellSpace / 2;
      const y = rows.indexOf(d.y) * cellSize + cellSpace / 2;
      const rectSize = cellSize - cellSpace;

      ctx.fillStyle = colorFn(d.value);
      ctx.fillRect(x, y, rectSize, rectSize);

      if (showPercentIdentities) {
        const textColor = tinycolor(colorFn(d.value)).isLight()
          ? "#000"
          : "#fff";
        ctx.fillStyle = textColor;
        ctx.textAlign = "center";
        ctx.textBaseline = "middle";
        ctx.font = `${annotation_font_size}px ${plotFontMonospace.family}`;
        ctx.fillText(
          d.value.toFixed(roundTo),
          x + rectSize / 2,
          y + rectSize / 2,
        );
      }
    }
    ctx.restore();

    // Titles
    if (showTitles) {
      ctx.fillStyle = "black";
      ctx.textAlign = "center";
      ctx.textBaseline = "top";
      ctx.font = `Bold 20px ${titleFont.family}`;
      ctx.fillText(title, width / 2, margin.top - 20);
      ctx.font = `20px ${titleFont.family}`;
      ctx.fillText(subtitle, width / 2, margin.top);
    }

    const axisGap = 5;

    if (axis_labels) {
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
        ctx.rotate((axlabel_xrotation * Math.PI) / 180);
        ctx.fillStyle = "black";
        ctx.textAlign = "right";
        ctx.textBaseline = "middle";
        ctx.font = `${Math.max(
          axlabel_xfontsize * transform.k,
          axlabel_xfontsize,
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
        ctx.rotate((axlabel_yrotation * Math.PI) / 180);
        ctx.fillStyle = "black";
        ctx.textAlign = "right";
        ctx.textBaseline = "middle";
        ctx.font = `${Math.max(
          axlabel_yfontsize * transform.k,
          axlabel_yfontsize,
        )}px ${plotFontMonospace.family}`;
        ctx.fillText(txt, 0, 0);
        ctx.restore();
      }
    }

    // Colorbar
    if (showscale) {
      const positionX = width - cbarWidth - margin.right;
      const gradient = ctx.createLinearGradient(
        minValue,
        margin.top + cbarHeight,
        minValue,
        margin.top,
      );

      for (const [stop, color] of colorScale) {
        gradient.addColorStop(stop, color);
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
  }, [
    transform,
    filteredData,
    colorFn,
    scale,
    tickValues,
    cellSize,
    cellSpace,
    showPercentIdentities,
    roundTo,
    showTitles,
    title,
    subtitle,
    annotation_font_size,
    axlabel_xfontsize,
    axlabel_yfontsize,
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
    minValue,
  ]);

  React.useEffect(() => {
    drawCanvas();
  }, [drawCanvas]);

  React.useEffect(() => {
    const canvas = canvasRef.current;
    if (!canvas) return;

    const zoom = d3
      .zoom()
      .scaleExtent([1, 5])
      .translateExtent([
        [-margin.left, -margin.top],
        [width, height],
      ])
      .on("zoom", (event) => setTransform(event.transform));

    d3.select(canvas).call(
      zoom as unknown as d3.ZoomBehavior<HTMLCanvasElement, unknown>,
    );

    return () => {
      d3.select(canvas).on(".zoom", null);
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

    const cell = filteredData.find((d) => d.x === dataX && d.y === dataY);

    if (cell) {
      setTooltipData({
        x,
        y,
        value: cell.value,
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
            <dt>Percent ID:</dt>
            <dd>{tooltipData.value.toFixed(2)}%</dd>
          </div>
        </dl>
      )}
    </div>
  );
};
