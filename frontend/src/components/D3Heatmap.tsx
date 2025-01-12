import * as d3 from "d3";
import React, { useEffect, useRef } from "react";
import tinycolor from "tinycolor2";
import type { ColorScaleArray } from "../colorScales";

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
  const [svgTransform, setSvgTransform] = React.useState({});

  useEffect(() => {
    if (!svgRef.current) return;
    const d3Svg = d3.select(svgRef.current as Element);
    d3Svg.selectAll("*").remove();

    const margin = { top: 60, right: 60, bottom: 60, left: 60 };
    // const margin = { top: 0, right: 0, bottom: 0, left: 0 };
    const w = width - margin.left - margin.right;
    const h = height - margin.top - margin.bottom;

    const n = tickText.length;

    const colorFn = createD3ColorScale(colorScale, minVal, maxVal);

    const g = d3
      .select(svgRef.current)
      .append("g")
      .attr("transform", `translate(${margin.left}, ${margin.top})`);

    const cellW = w / n;
    const cellH = h / n;
    const cellOffset = cellSpace > 0 ? cellSpace / 2 : 0;
    const labelOffset = 1;
    const fontSizeMin = 1;
    const fontSizeMax = 16;
    const fontSizeFactor = 0.01;
    const fontSize =
      fontSizeMin +
      (fontSizeMax - fontSizeMin) /
        (fontSizeMin + fontSizeFactor * data.length);

    const groups = g
      .selectAll("g")
      .data(data.filter((d) => Number(d.value)))
      .join("g")
      .attr("transform", (d) => `translate(${d.x * cellW}, ${d.y * cellH})`);

    groups
      .append("rect")
      .attr("width", cellW - cellSpace)
      .attr("height", cellH - cellSpace)
      .attr("x", cellOffset)
      .attr("y", cellOffset)
      .attr("fill", (d) => colorFn(d.value));

    groups
      .append("text")
      .attr("x", cellW / 2)
      .attr("y", cellH / 2)
      .attr("dy", ".35em")
      .attr("text-anchor", "middle")
      .attr("font-family", "Roboto Mono")
      .attr("font-size", `${fontSize}px`)
      .text((d) => d.value)
      .attr("fill", (d) =>
        tinycolor(colorFn(d.value)).isLight() ? "#000" : "#fff",
      );

    d3Svg.call(
      d3
        .zoom()
        .extent([
          [0, 0],
          [width, height],
        ])
        .scaleExtent([1, 25])
        .translateExtent([
          [0 - margin.left, 0 - margin.top],
          [width, height],
        ])
        .on("zoom", ({ transform }: d3.D3ZoomEvent<SVGSVGElement, unknown>) => {
          setSvgTransform(transform as object);
          g.attr("transform", transform.toString());
        }),
    );

    // Simple x-axis labels
    const xAxisG = g.append("g").attr("transform", `translate(0, ${h})`);
    xAxisG
      .selectAll("text")
      .data(tickText)
      .join("text")
      .attr("y", cellH)
      .attr("text-anchor", "middle")
      .attr("font-family", "Roboto Mono")
      .attr("font-size", `${fontSize}px`)
      .attr("font-weight", "bold")
      .attr("text-anchor", "end")
      .attr(
        "transform",
        (_, i) =>
          `translate(${i * cellW - cellW / 2}, ${labelOffset}) rotate(-90)`,
      )
      .text((txt) => txt);

    // Simple y-axis labels
    const yAxisG = g.append("g");
    yAxisG
      .selectAll("text")
      .data(tickText)
      .join("text")
      .attr("x", -labelOffset)
      .attr("y", (_, i) => i * cellH + cellH / 2)
      .attr("dominant-baseline", "middle")
      .attr("text-anchor", "end")
      .attr("font-family", "Roboto Mono")
      .attr("font-size", `${fontSize}px`)
      .attr("font-weight", "bold")
      .text((txt) => txt);
  }, [data, tickText, colorScale, minVal, maxVal, width, height, cellSpace]);

  console.log(svgTransform);

  return (
    <svg
      style={{ background: "#fff" }}
      ref={svgRef}
      width={"100%"}
      height={"100%"}
    />
  );
};
