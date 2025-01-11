import * as d3 from "d3";
import { useEffect, useRef } from "react";
import { type ColorScaleArray, colorScales } from "../colorScales";

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

  // Build a piecewise linear scale
  return d3
    .scaleLinear<string>()
    .domain(domain)
    .range(range)
    .interpolate(d3.interpolateRgb);
}

export const D3Heatmap = ({
  data,
  tickText,
  colorScale,
  minVal = 0,
  maxVal = 100,
  width = 500,
  height = 500,
  cellSpace,
}: D3HeatmapProps) => {
  const svgRef = useRef<SVGSVGElement | null>(null);

  useEffect(() => {
    if (!svgRef.current) return;
    d3.select(svgRef.current).selectAll("*").remove();

    const margin = { top: 60, right: 60, bottom: 60, left: 60 };
    const w = width - margin.left - margin.right;
    const h = height - margin.top - margin.bottom;

    // Number of rows/cols from tickText
    const n = tickText.length;

    const colorFn = createD3ColorScale(colorScale, minVal, maxVal);

    const g = d3
      .select(svgRef.current)
      .append("g")
      .attr("transform", `translate(${margin.left}, ${margin.top})`);

    const cellW = w / n;
    const cellH = h / n;

    // Draw each cell
    g.selectAll("rect")
      .data(data.filter((d) => Number(d.value)))
      .join("rect")
      .attr("x", (d) => d.x * cellW + cellSpace / 2)
      .attr("y", (d) => d.y * cellH + cellSpace / 2)
      .attr("width", cellW - cellSpace)
      .attr("height", cellH - cellSpace)
      .attr("fill", (d) => colorFn(d.value))
      .append("title")
      .text((d) => `Value: ${d.value.toFixed(1)}%`);

    // Simple x-axis labels
    const xAxisG = g.append("g").attr("transform", `translate(0, ${h})`);
    xAxisG
      .selectAll("text")
      .data(tickText)
      .join("text")
      .attr("x", (_, i) => i * cellW + cellW / 2)
      .attr("y", 20)
      .attr("text-anchor", "middle")
      .text((txt) => txt);

    // Simple y-axis labels
    const yAxisG = g.append("g");
    yAxisG
      .selectAll("text")
      .data(tickText)
      .join("text")
      .attr("x", -5)
      .attr("y", (_, i) => i * cellH + cellH / 2)
      .attr("dominant-baseline", "middle")
      .attr("text-anchor", "end")
      .text((txt) => txt);
  }, [data, tickText, colorScale, minVal, maxVal, width, height, cellSpace]);

  return <svg ref={svgRef} width={width} height={height} />;
};
