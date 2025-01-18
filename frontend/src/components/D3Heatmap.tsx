import * as d3 from "d3";
import React from "react";
import tinycolor from "tinycolor2";
import type { ColorScaleArray } from "../colorScales";
import { plotFontMonospace, plotFontSansSerif } from "../constants";
import { useHeatmapRef } from "../hooks/useHeatmapRef";
import type { HeatmapSettings } from "../plotTypes";
import { ColorLegend } from "./ColorLegend";

interface HeatmapCell {
  x: number;
  y: number;
  value: number;
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
  // axis_labels = true,
  cellSpace,
  roundTo,
  showPercentIdentities = true,
  showscale,
  cbarWidth,
  cbarHeight,
  axlabel_xfontsize,
  axlabel_yfontsize,
  axlabel_xrotation = 90,
  axlabel_yrotation = 180,
  tempHeatmapComponent,
  titleFont,
  showTitles = true,
  title = "",
  subtitle = "",
}: {
  data: HeatmapCell[];
  tickText: string[];
  colorScale: ColorScaleArray;
  minVal?: number;
  maxVal?: number;
  width?: number;
  height?: number;
  cellSpace: number;
  roundTo: number;
  showPercentIdentities: boolean;
  showscale?: boolean;
  cbarWidth: number;
  cbarHeight: number;
  // axis_labels: boolean;
  axlabel_xfontsize: number;
  axlabel_yfontsize: number;
  axlabel_xrotation: number;
  axlabel_yrotation: number;
  tempHeatmapComponent: "canvas" | "svg" | "plotly";
  titleFont: HeatmapSettings["titleFont"];
  showTitles: boolean;
  title: string;
  subtitle: string;
}) => {
  const svgRef = useHeatmapRef() as React.MutableRefObject<SVGSVGElement>;
  const [svgTransform, setSvgTransform] = React.useState({});
  console.log(svgTransform);

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
    if (!svgRef.current) return;
    const d3Svg = d3.select(svgRef.current as Element);
    d3Svg.selectAll("*").remove();

    const plotFont =
      titleFont === "Monospace" ? plotFontMonospace : plotFontSansSerif;

    const size = Math.min(width, height);
    const margin = { top: 60, right: 60, bottom: 60, left: 60 };
    const w = size - margin.left - margin.right;
    const h = size - margin.top - margin.bottom;

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

    if (showPercentIdentities) {
      groups
        .append("text")
        .attr("x", cellW / 2)
        .attr("y", cellH / 2)
        .attr("dy", ".35em")
        .attr("text-anchor", "middle")
        .attr("font-family", "Roboto Mono")
        .attr("font-size", `${fontSize}px`)
        .text((d) => d.value.toFixed(roundTo))
        .attr("fill", (d) =>
          tinycolor(colorFn(d.value)).isLight() ? "#000" : "#fff",
        );
    }

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
    // titles
    if (showTitles) {
      console.log("Title props:", {
        showTitles,
        title,
        width,
        height,
        fontSize,
        plotFont,
      });
      g.append("text")
        .attr("text-anchor", "middle") // try middle instead of end
        .attr("font-family", plotFont.family)
        .attr("font-size", `${fontSize}px`)
        .attr("x", width / 2)
        .attr("y", margin.top - 20) // try margin.top instead of height
        .text(`${title}`);

      g.append("text")
        .attr("text-anchor", "middle") // try middle instead of end
        .attr("font-family", plotFont.family)
        .attr("font-size", `${fontSize}px`)
        .attr("x", width / 2)
        .attr("y", margin.top - labelOffset * 2) // try margin.top instead of height
        .text(`${subtitle}`);
    }
    // x-axis labels
    g.append("g")
      .attr("transform", `translate(0, ${h})`)
      .selectAll("text")
      .data(tickText)
      .join("text")
      .attr("x", (_, i) => i * cellW + cellW / 2)
      .attr("y", labelOffset)
      .attr("dy", "1em")
      .attr("text-anchor", "end")
      .attr("font-family", plotFont.family)
      .attr("font-size", axlabel_xfontsize)
      .text((txt) => txt)
      .attr(
        "transform",
        (_, i) =>
          `rotate(${axlabel_xrotation}, ${i * cellW + cellW / 2}, ${labelOffset})`,
      );

    // y-axis labels
    g.append("g")
      .selectAll("text")
      .data(tickText)
      .join("text")
      .attr("x", -labelOffset)
      .attr("y", (_, i) => i * cellH + cellH / 2)
      .attr("dominant-baseline", "middle")
      .attr("text-anchor", "end")
      .attr("font-family", plotFont.family)
      .attr("font-size", axlabel_yfontsize)
      .text((txt) => txt)
      .attr(
        "transform",
        (_, i) =>
          `rotate(${axlabel_yrotation}, ${-labelOffset}, ${i * cellH + cellH / 2})`,
      );

    const scaleBarX = width - cbarWidth - margin.left - margin.right;

    const gradientGroup = g
      .append("g")
      .attr("transform", `translate(${scaleBarX},${margin.top})`);
    const defs = d3Svg.append("defs");
    const gradientId = "heatmap-svg-linear-gradient";
    const gradient = defs
      .append("linearGradient")
      .attr("id", gradientId)
      .attr("x1", "0%")
      .attr("y1", "100%")
      .attr("x2", "0%")
      .attr("y2", "0%");

    for (const { stop, color } of gradientStops) {
      gradient
        .append("stop")
        .attr("offset", `${stop * 100}%`)
        .attr("stop-color", color);
    }

    // colorScale.map(([stop, color]) => {
    //   gradient
    //     .append("stop")
    //     .attr("offset", `${stop * 100}%`)
    //     .attr("stop-color", color);
    // });

    // const legendHeight = cbarHeight - margin.top - margin.bottom;
    // const legendWidth = cbarWidth - margin.left - margin.right;

    gradientGroup
      .append("rect")
      .attr("width", cbarWidth)
      .attr("height", cbarHeight)
      .attr("fill", `url(#${gradientId})`);

    const axis = d3.axisRight(scale).tickValues(tickValues).tickSize(6);

    g.append("g")
      .attr("transform", `translate(${scaleBarX + cbarWidth},${margin.top})`)
      .call(axis)
      .call((g) => g.select(".domain").remove())
      .call((g) => g.selectAll(".tick line").attr("stroke-opacity", 0.5))
      .call((g) =>
        g
          .selectAll(".tick text")
          .attr("font-size", "10px")
          .attr("font-family", "Roboto Mono")
          .attr("dx", "0.5em"),
      );
  }, [
    tickValues,
    scale,
    gradientStops,
    cbarHeight,
    cbarWidth,
    svgRef.current,
    data,
    tickText,
    colorScale,
    minVal,
    maxVal,
    width,
    height,
    cellSpace,
    roundTo,
    showPercentIdentities,
    axlabel_xfontsize,
    axlabel_yfontsize,
    axlabel_xrotation,
    axlabel_yrotation,
    titleFont,
    showTitles,
    title,
    subtitle,
  ]);

  return (
    <svg
      xmlns="http://www.w3.org/2000/svg"
      xmlnsXlink="http://www.w3.org/1999/xlink"
      style={{ background: "#fff" }}
      ref={svgRef}
      width={"100%"}
      height={"100%"}
    />
  );
};
