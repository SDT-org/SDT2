import * as d3 from "d3";
import React from "react";
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
  showscale?: boolean;
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
  // showscale = true,
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

  const filteredData = data.filter((d) => Number(d.value));

  const gradientStops = React.useMemo(
    () => colorScale.map(([stop, color]) => ({ stop, color })),
    [colorScale],
  );
  const scale = React.useMemo(
    () => d3.scaleLinear().domain([maxVal, minVal]).range([0, cbarHeight]),
    [minVal, maxVal, cbarHeight],
  );
  const ticks = 5;
  const tickValues = React.useMemo(() => scale.ticks(ticks), [scale]);

  React.useEffect(() => {
    const canvas = canvasRef.current;
    const plotFont =
      titleFont === "Monospace" ? plotFontMonospace : plotFontSansSerif;
    if (!canvas) return;

    const ctx = canvas.getContext("2d");
    if (!ctx) return;
    ctx.setTransform(transform.k, 0, 0, transform.k, transform.x, transform.y);
    ctx.clearRect(0, 0, canvas.width, canvas.height);
    ctx.save();
    // Calculate  dimensions with forced aspect of same sidedness
    const size = Math.min(width, height);
    const pixelRatio = window.devicePixelRatio || 1;

    // Set canvas to be square
    canvas.width = size * pixelRatio;
    canvas.height = size * pixelRatio;
    canvas.style.width = `${size}px`;
    canvas.style.height = `${size}px`;
    ctx.scale(pixelRatio, pixelRatio);

    // Scale margins with zoom to match hover detection
    const margin = {
      top: 60 * transform.k,
      right: 60 * transform.k,
      bottom: 60 * transform.k,
      left: 60 * transform.k,
    };
    // set plot size
    const plotSize = size - margin.left - margin.right;
    const n = tickText.length;
    const cellSize = plotSize / n;

    const colorFn = createD3ColorScale(colorScale, minVal, maxVal);
    const rows = [...new Set(filteredData.map((d) => d.x))];
    const cols = [...new Set(filteredData.map((d) => d.y))];

    ctx.clearRect(0, 0, size, size);
    ctx.save();
    //un does panning
    ctx.translate(transform.x, transform.y);
    ctx.scale(transform.k, transform.k);

    //descale the
    ctx.translate(margin.left / transform.k, margin.top / transform.k);

    // Draw cells
    for (const d of filteredData) {
      const x = cols.indexOf(d.x) * cellSize + cellSpace / 2;
      const y = rows.indexOf(d.y) * cellSize + cellSpace / 2;
      const rectSize = cellSize - cellSpace;

      // draw rects fill color
      ctx.fillStyle = colorFn(d.value);
      ctx.fillRect(x, y, rectSize, rectSize);

      ctx.save();
      const fontSizeMin = 1;
      const fontSizeMax = 16;
      const fontSizeFactor = 0.01;
      const fontSize = Math.min(
        fontSizeMax,
        Math.max(fontSizeMin, rectSize + fontSizeFactor * data.length),
      );

      ctx.fillStyle = "black";
      ctx.textAlign = "center";
      ctx.textBaseline = "middle";
      ctx.font = `${fontSize}px ${plotFont.family}`;
      if (showPercentIdentities) {
        ctx.fillText(
          `${d.value.toFixed(roundTo)}`,
          x + rectSize / 2,
          y + rectSize / 2,
        );
      }
      ctx.restore();
    }

    // titles
    ctx.fillStyle = "black";
    ctx.textAlign = "center";
    ctx.textBaseline = "top";
    if (showTitles) {
      ctx.font = `${axlabel_xfontsize}px ${plotFont.family}`;
      ctx.fillText(`${title}`, plotSize / 2, -margin.top + 20);
      ctx.fillText(`${subtitle}`, plotSize / 2, -margin.top + 40);
    }

    // x axis labels
    ctx.fillStyle = "black";
    ctx.textBaseline = "middle";
    ctx.font = `${axlabel_xfontsize}px ${plotFont.family}`;
    ctx.textAlign = "right";
    for (let i = 0; i < tickText.length; i++) {
      const txt = tickText[i];
      if (txt === undefined) continue;
      ctx.save();
      ctx.translate(i * cellSize + cellSize / 2, plotSize + 20);
      ctx.rotate((axlabel_xrotation * Math.PI) / 180);
      ctx.fillText(txt, 0, 0);
      ctx.restore();
    }

    // y axis labels
    ctx.fillStyle = "black";
    ctx.textBaseline = "middle";
    ctx.font = `${axlabel_yfontsize}px ${plotFont.family}`;
    ctx.textAlign = "right";
    for (let i = 0; i < tickText.length; i++) {
      const txt = tickText[i];
      if (txt === undefined) continue;
      ctx.save();
      ctx.translate(-20, i * cellSize + cellSize / 2);
      ctx.rotate((axlabel_yrotation * Math.PI) / 180);
      ctx.fillText(txt, 0, 0);
      ctx.restore();
    }

    // colorbar
    const positionX = width - cbarWidth - margin.left - margin.right;

    const gradient = ctx.createLinearGradient(
      0,
      cbarHeight + margin.top,
      0,
      margin.top,
    );

    for (const { stop, color } of gradientStops) {
      gradient.addColorStop(stop, color);
    }

    ctx.fillStyle = gradient;
    ctx.fillRect(positionX, margin.top, cbarWidth, cbarHeight);

    const adjustedScale = scale
      .copy()
      .range([margin.top, cbarHeight + margin.top]);

    ctx.fillStyle = "#000";
    ctx.font = "10px 'Roboto Mono'";
    ctx.textAlign = "left";
    ctx.textBaseline = "middle";

    for (const tick of tickValues) {
      const y = adjustedScale(tick);
      ctx.fillText(tick.toFixed(2), positionX + cbarWidth + 5, y);
      ctx.beginPath();
      ctx.moveTo(positionX + cbarWidth, y);
      ctx.lineTo(positionX + cbarWidth + 6, y);
      ctx.strokeStyle = "#000";
      ctx.stroke();
    }
  }, [
    cbarHeight,
    cbarWidth,
    scale,
    gradientStops,
    tickValues,
    canvasRef.current,
    data,
    tickText,
    colorScale,
    minVal,
    maxVal,
    width,
    height,
    cellSpace,
    roundTo,
    transform,
    axlabel_xrotation,
    axlabel_yrotation,
    axlabel_xfontsize,
    axlabel_yfontsize,
    titleFont,
    showPercentIdentities,
    showTitles,
    title,
    subtitle,
    filteredData,
  ]);

  React.useEffect(() => {
    const canvas = canvasRef.current;
    if (!canvas) return;

    const zoom = d3
      .zoom()
      .scaleExtent([0.5, 5])
      .on("zoom", (event) => {
        setTransform(event.transform);
      });

    // Type assertion for d3.zoom() call
    d3.select(canvas).call(
      zoom as unknown as (
        selection: d3.Selection<HTMLCanvasElement, unknown, null, undefined>,
      ) => void,
    );
  }, [canvasRef]);

  const handleMouseMove = (event: React.MouseEvent<HTMLCanvasElement>) => {
    const canvas = canvasRef.current;
    if (!canvas) return;
    const size = Math.min(width, height);

    const rect = canvas.getBoundingClientRect();
    const x = event.clientX - rect.left;
    const y = event.clientY - rect.top;

    // Keep your correct margin scaling
    const margin = {
      top: 60 * transform.k,
      right: 60 * transform.k,
      bottom: 60 * transform.k,
      left: 60 * transform.k,
    };

    // Match drawing calcs but with scaled margins
    const plotSize = size - margin.left - margin.right;
    const n = tickText.length;
    const cellSize = plotSize / n;

    //calcualte the data x and y by subtracting pan and zoom
    const dataX = Math.floor(
      (x - margin.left - transform.x) / (cellSize * transform.k),
    );
    const dataY = Math.floor(
      (y - margin.top - transform.y) / (cellSize * transform.k),
    );

    const cell = filteredData.find((d) => d.x === dataX && d.y === dataY);

    if (cell) {
      setTooltipData({ x, y, value: cell.value });
    } else {
      setTooltipData(null);
    }
  };

  const handleMouseLeave = () => {
    setTooltipData(null);
  };

  return (
    <div style={{ position: "relative" }}>
      <canvas
        ref={canvasRef}
        style={{ background: "#fff" }}
        width={width}
        height={height}
        onMouseMove={handleMouseMove}
        onMouseLeave={handleMouseLeave}
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
