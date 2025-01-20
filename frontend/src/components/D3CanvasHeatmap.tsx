import * as d3 from "d3";
import React from "react";
import tinycolor from "tinycolor2";
import type { ColorScaleArray } from "../colorScales";
import { plotFontMonospace, plotFontSansSerif } from "../constants";
import { useHeatmapRef } from "../hooks/useHeatmapRef";
import type { HeatmapSettings } from "../plotTypes";

interface HeatmapCell {
  x: number;
  y: number;
  value: number;
}

interface D3HeatmapProps {
  data: HeatmapCell[];
  tickText: string[];
  colorScale: ColorScaleArray;
  minVal?: number;
  maxVal?: number;
  width?: number;
  height?: number;
  cellSpace: number;
  roundTo: number;
  cbarWidth?: number;
  cbarHeight?: number;
  axlabel_xfontsize: number;
  axlabel_yfontsize: number;
  axlabel_xrotation: number;
  axlabel_yrotation: number;
  titleFont: HeatmapSettings["titleFont"];
  showPercentIdentities: boolean;
  showTitles: boolean;
  title: string;
  subtitle: string;
  showscale?: boolean;
}

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

const defaultMargin = {
  top: 60,
  right: 60,
  bottom: 60,
  left: 60,
} as const;

export const D3CanvasHeatmap = ({
  data,
  tickText,
  colorScale,
  minVal = 0,
  maxVal = 100,
  width = 500,
  height = 500,
  cellSpace,
  roundTo,
  cbarHeight = 200,
  cbarWidth = 60,
  axlabel_xfontsize = 15,
  axlabel_yfontsize = 15,
  axlabel_xrotation = 0,
  axlabel_yrotation = 0,
  titleFont,
  showPercentIdentities = true,
  showTitles = true,
  title = "DOODOO",
  subtitle = "",
}: D3HeatmapProps) => {
  const canvasRef =
    useHeatmapRef() as React.MutableRefObject<HTMLCanvasElement>;
  const [transform, setTransform] = React.useState(d3.zoomIdentity);
  const [tooltipData, setTooltipData] = React.useState<{
    x: number;
    y: number;
    value: number;
  } | null>(null);

  const filteredData = data.filter((d: HeatmapCell) => Number(d.value));
  const size = Math.min(width, height);
  const plotSize = size - defaultMargin.left - defaultMargin.right;
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

  const drawCanvas = React.useCallback(() => {
    const canvas = canvasRef.current;
    if (!canvas) return;

    const ctx = canvas.getContext("2d");
    if (!ctx) return;

    const pixelRatio = window.devicePixelRatio || 1;
    canvas.width = size * pixelRatio;
    canvas.height = size * pixelRatio;
    canvas.style.width = `${size}px`;
    canvas.style.height = `${size}px`;

    ctx.clearRect(0, 0, canvas.width, canvas.height);
    ctx.scale(pixelRatio, pixelRatio);

    // set zoom transform
    ctx.save();
    ctx.translate(
      transform.x + defaultMargin.left,
      transform.y + defaultMargin.top,
    );
    ctx.scale(transform.k, transform.k);
    //indexz data
    const rows = [...new Set(filteredData.map((d: HeatmapCell) => d.x))];
    const cols = [...new Set(filteredData.map((d: HeatmapCell) => d.y))];

    // Draw cells
    for (const d of filteredData) {
      const x = cols.indexOf(d.x) * cellSize + cellSpace / 2;
      const y = rows.indexOf(d.y) * cellSize + cellSpace / 2;
      const rectSize = cellSize - cellSpace;

      ctx.fillStyle = colorFn(d.value);
      ctx.fillRect(x, y, rectSize, rectSize);
      // Percent ids
      if (showPercentIdentities) {
        const fontSize = Math.min(
          16,
          Math.max(1, rectSize + 0.01 * data.length),
        );
        const textColor = tinycolor(colorFn(d.value)).isLight()
          ? "#000"
          : "#fff";
        ctx.fillStyle = textColor;
        ctx.textAlign = "center";
        ctx.textBaseline = "middle";
        ctx.font = `${fontSize}px ${titleFont === "Monospace" ? plotFontMonospace.family : plotFontSansSerif.family}`;
        ctx.fillText(
          d.value.toFixed(roundTo),
          x + rectSize / 2,
          y + rectSize / 2,
        );
      }
    }
    ctx.restore();

    // font settings
    const plotFont =
      titleFont === "Monospace"
        ? plotFontMonospace.family
        : plotFontSansSerif.family;

    // Titles
    if (showTitles) {
      ctx.fillStyle = "black";
      ctx.textAlign = "center";
      ctx.textBaseline = "top";
      ctx.font = `${axlabel_xfontsize}px ${plotFont}`;
      ctx.fillText(title, size / 2, 20);
      ctx.fillText(subtitle, size / 2, 40);
    }

    // X-axis labels
    for (const [i, txt] of tickText.entries()) {
      if (txt === undefined) continue;
      ctx.save();
      ctx.translate(
        defaultMargin.left + //margin offset
          i * cellSize * transform.k + //current tick position
          (cellSize * transform.k) / 2 + //
          transform.x,
        size - defaultMargin.bottom + 20,
      );
      ctx.rotate((axlabel_xrotation * Math.PI) / 180);
      ctx.fillStyle = "black";
      ctx.textAlign = "right";
      ctx.textBaseline = "middle";
      ctx.font = `${axlabel_xfontsize}px ${plotFont}`;
      ctx.fillText(txt, 0, 0);
      ctx.restore();
    }

    // Y-axis labels
    for (const [i, txt] of tickText.entries()) {
      if (txt === undefined) continue;
      ctx.save();
      ctx.translate(
        defaultMargin.left - 20, // X position: margin offset with spacing for labels
        defaultMargin.top + // Y start position (top margin)
          i * cellSize * transform.k + // Y position for current tick (accounting for zoom)
          (cellSize * transform.k) / 2 + // Center vertically within cell
          transform.y, // Y pan offset
      );
      ctx.rotate((axlabel_yrotation * Math.PI) / 180);
      ctx.fillStyle = "black";
      ctx.textAlign = "right";
      ctx.textBaseline = "middle";
      ctx.font = `${axlabel_yfontsize}px ${plotFont}`;
      ctx.fillText(txt, 0, 0);
      ctx.restore();
    }

    // Colorbar
    const positionX = size - cbarWidth - defaultMargin.right;
    const gradient = ctx.createLinearGradient(
      0,
      defaultMargin.top + cbarHeight,
      0,
      defaultMargin.top,
    );

    for (const [stop, color] of colorScale) {
      gradient.addColorStop(stop, color);
    }

    ctx.fillStyle = gradient;
    ctx.fillRect(positionX, defaultMargin.top, cbarWidth, cbarHeight);

    ctx.fillStyle = "black";
    ctx.textAlign = "left";
    ctx.textBaseline = "middle";
    ctx.font = "10px 'Roboto Mono'";

    for (const tick of tickValues) {
      const y = scale(tick) + defaultMargin.top;
      ctx.fillText(tick.toFixed(2), positionX + cbarWidth + 5, y);
      ctx.beginPath();
      ctx.moveTo(positionX + cbarWidth, y);
      ctx.lineTo(positionX + cbarWidth + 6, y);
      ctx.stroke();
    }
  }, [
    transform,
    filteredData,
    colorFn,
    scale,
    tickValues,
    size,
    cellSize,
    cellSpace,
    showPercentIdentities,
    roundTo,
    showTitles,
    title,
    subtitle,
    axlabel_xfontsize,
    axlabel_yfontsize,
    axlabel_xrotation,
    axlabel_yrotation,
    titleFont,
    tickText,
    cbarWidth,
    cbarHeight,
    colorScale,
    data.length,
    canvasRef,
  ]);

  React.useEffect(() => {
    drawCanvas();
  }, [drawCanvas]);

  React.useEffect(() => {
    const canvas = canvasRef.current;
    if (!canvas) return;

    const zoom = d3
      .zoom()
      .scaleExtent([0.5, 5])
      .on("zoom", (event) => {
        setTransform(event.transform);
      });

    d3.select(canvas).call(
      zoom as unknown as d3.ZoomBehavior<HTMLCanvasElement, unknown>,
    );

    return () => {
      d3.select(canvas).on(".zoom", null);
    };
  }, [canvasRef]);

  const handleMouseMove = (event: React.MouseEvent<HTMLCanvasElement>) => {
    const canvas = canvasRef.current;
    if (!canvas) return;

    const rect = canvas.getBoundingClientRect();
    const x = event.clientX - rect.left;
    const y = event.clientY - rect.top;

    const dataX = Math.floor(
      (x - defaultMargin.left - transform.x) / (cellSize * transform.k),
    );
    const dataY = Math.floor(
      (y - defaultMargin.top - transform.y) / (cellSize * transform.k),
    );

    const cell = filteredData.find(
      (d: HeatmapCell) => d.x === dataX && d.y === dataY,
    );

    if (cell) {
      setTooltipData({ x, y, value: cell.value });
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
        <div
          style={{
            position: "absolute",
            left: tooltipData.x,
            top: tooltipData.y,
            background: "white",
            border: "1px solid black",
            padding: "20px",
            pointerEvents: "none",
          }}
        >
          {tooltipData.value.toFixed(2)}
        </div>
      )}
    </div>
  );
};
