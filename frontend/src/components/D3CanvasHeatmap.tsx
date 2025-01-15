import * as d3 from "d3";
import { useEffect, useRef, useState } from "react";
import type { ColorScaleArray } from "../colorScales";
import { plotFontMonospace, plotFontSansSerif } from "../constants";
import type { HeatmapSettings } from "../plotTypes";
import { ColorLegend } from "./ColorLegend";

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
  showscale?: boolean;
  cbarWidth?: number;
  cbarHeight?: number;
  axlabel_xfontsize: number;
  axlabel_yfontsize: number;
  axlabel_xrotation: number;
  axlabel_yrotation: number;
  titleFont: HeatmapSettings["titleFont"];
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
  showscale = true,
  cbarHeight = 200,
  cbarWidth = 60,
  axlabel_xfontsize = 15,
  axlabel_yfontsize = 15,
  axlabel_xrotation = 0,
  axlabel_yrotation = 0,
  titleFont,
}: D3HeatmapProps) => {
  const canvasRef = useRef<HTMLCanvasElement>(null);
  const [transform, setTransform] = useState(d3.zoomIdentity);
  const [tooltipData, setTooltipData] = useState<{
    x: number;
    y: number;
    value: number;
  } | null>(null);

  useEffect(() => {
    const canvas = canvasRef.current;
    const plotFont =
      titleFont === "Monospace" ? plotFontMonospace : plotFontSansSerif;
    if (!canvas) return;

    const ctx = canvas.getContext("2d");
    if (!ctx) return;

    const pixelRatio = window.devicePixelRatio || 1;
    canvas.width = width * pixelRatio;
    canvas.height = height * pixelRatio;
    canvas.style.width = `${width}px`;
    canvas.style.height = `${height}px`;
    ctx.scale(pixelRatio, pixelRatio);

    const margin = { top: 60, right: 60, bottom: 60, left: 60 };
    const w = width - margin.left - margin.right;
    const h = height - margin.top - margin.bottom;

    const n = tickText.length;
    const colorFn = createD3ColorScale(colorScale, minVal, maxVal);

    const cellW = w / n;
    const cellH = h / n;

    const rows = [...new Set(data.map((d) => d.x))];
    const cols = [...new Set(data.map((d) => d.y))];

    ctx.clearRect(0, 0, width, height);
    ctx.save();

    ctx.translate(transform.x, transform.y);
    ctx.scale(transform.k, transform.k);
    ctx.translate(margin.left, margin.top);

    const filteredData = data.filter((d) => Number(d.value));
    for (const d of filteredData) {
      const x = cols.indexOf(d.x) * cellW + cellSpace / 2;
      const y = rows.indexOf(d.y) * cellH + cellSpace / 2;
      const rectW = cellW - cellSpace;
      const rectH = cellH - cellSpace;

      ctx.fillStyle = colorFn(d.value);
      ctx.fillRect(x, y, rectW, rectH);
      ctx.save();
      ctx.fillStyle = "black";
      ctx.textAlign = "center";
      ctx.textBaseline = "middle";

      ctx.fillText(`${d.value}`, x + rectW / 2, y + rectH / 2);

      ctx.restore();
    }

    ctx.fillStyle = "black";
    ctx.textAlign = "center";
    ctx.textBaseline = "middle";
    ctx.font = `${axlabel_xfontsize}px ${plotFont.family}`;
    ctx.textAlign = "right";
    for (let i = 0; i < tickText.length; i++) {
      const txt = tickText[i];
      if (txt === undefined) continue;
      ctx.save();
      ctx.translate(i * cellW + cellW / 2, h + 20);
      ctx.rotate((axlabel_xrotation * Math.PI) / 180);
      ctx.fillText(txt, 0, 0);
      ctx.restore();
    }
    ctx.font = `${axlabel_yfontsize}px ${plotFont.family}`;
    ctx.textAlign = "right";
    for (let i = 0; i < tickText.length; i++) {
      const txt = tickText[i];
      if (txt === undefined) continue;
      ctx.save();
      ctx.translate(-5, i * cellH + cellH / 2);
      ctx.rotate((axlabel_yrotation * Math.PI) / 180);
      ctx.fillText(txt, 0, 0);
      ctx.restore();
    }
  }, [
    data,
    tickText,
    colorScale,
    minVal,
    maxVal,
    width,
    height,
    cellSpace,
    transform,
    axlabel_xrotation,
    axlabel_yrotation,
    axlabel_xfontsize,
    axlabel_yfontsize,
    titleFont,
  ]);

  useEffect(() => {
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
  }, []);

  const handleMouseMove = (event: React.MouseEvent<HTMLCanvasElement>) => {
    const canvas = canvasRef.current;
    if (!canvas) return;

    const rect = canvas.getBoundingClientRect();
    const x = event.clientX - rect.left;
    const y = event.clientY - rect.top;

    const margin = { top: 60, right: 60, bottom: 60, left: 60 };
    const w = width - margin.left - margin.right;
    const h = height - margin.top - margin.bottom;
    const n = tickText.length;

    const cellW = w / n;
    const cellH = h / n;

    const dataX = Math.floor(
      (x - margin.left - transform.x) / (cellW * transform.k),
    );
    const dataY = Math.floor(
      (y - margin.top - transform.y) / (cellH * transform.k),
    );

    const cell = data.find((d) => d.x === dataX && d.y === dataY);

    if (cell) {
      setTooltipData({
        x: x + 10,
        y: y + 10,
        value: cell.value,
      });
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
      {showscale && (
        <ColorLegend
          colorScale={colorScale}
          minVal={minVal}
          maxVal={maxVal}
          position={{ x: width - cbarWidth - 20, y: height / 100 }}
          cbarHeight={cbarHeight}
          cbarWidth={cbarWidth}
          tempHeatmapComponent="canvas"
        />
      )}
      {tooltipData && (
        <div
          style={{
            position: "absolute",
            left: tooltipData.x,
            top: tooltipData.y,
            background: "white",
            border: "1px solid black",
            padding: "5px",
            pointerEvents: "none",
          }}
        >
          {tooltipData.value.toFixed(2)}
        </div>
      )}
    </div>
  );
};
