import * as d3 from "d3";
import React from "react";
import { plotFontMonospace } from "../constants";
import { getCellMetrics } from "../heatmapUtils";
import { useHeatmapRef } from "../hooks/useHeatmapRef";
import type { HeatmapRenderProps } from "./Heatmap";

export const D3Heatmap = ({
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
  cbarWidth,
  cbarHeight,
  axlabel_fontsize,
  axlabel_xrotation,
  axlabel_yrotation,
  titleFont,
  showTitles,
  title,
  showscale,
  axis_labels,
  margin,
}: HeatmapRenderProps) => {
  const svgRef = useHeatmapRef() as React.MutableRefObject<SVGSVGElement>;
  const [_, setSvgTransform] = React.useState({});

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
    if (!(svgRef.current && height && width)) return;
    const d3Svg = d3.select(svgRef.current as Element);
    d3Svg.selectAll("*").remove();

    const plotSize = Math.min(width, height);
    const plotWidth = plotSize - margin.left - margin.right;
    const plotHeight = plotSize - margin.top - margin.bottom;
    const cellSize = plotWidth / tickText.length;
    const cellMetrics = getCellMetrics(cellSize, cellSpace, roundTo + 3);

    const g = d3
      .select(svgRef.current)
      .append("g")
      .attr("transform", `translate(${margin.left}, ${margin.top})`);

    const axisGap = 5;

    const groups = g
      .selectAll("g")
      .data(data.filter((d) => Number(d.value)))
      .join("g")
      .attr(
        "transform",
        (d) => `translate(${d.x * cellSize}, ${d.y * cellSize})`,
      );

    groups
      .append("rect")
      .attr("width", Math.max(cellMetrics.cellSize, 1))
      .attr("height", Math.max(cellMetrics.cellSize, 1))
      .attr("x", cellMetrics.cellOffset)
      .attr("y", cellMetrics.cellOffset)
      .attr("fill", (d) => d.backgroundColor);

    if (showPercentIdentities) {
      groups
        .append("text")
        .attr("x", cellSize / 2)
        .attr("y", cellSize / 2)
        .attr("dy", ".35em")
        .attr("text-anchor", "middle")
        .attr("font-family", "Roboto Mono")
        .attr("font-size", `${cellMetrics.fontSize}px`)
        .text((d) => (d.value === 100 ? "100" : d.value.toFixed(roundTo)))
        .attr("fill", (d) => d.foregroundColor);
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

    if (showTitles) {
      g.append("text")
        .attr("text-anchor", "middle")
        .attr("font-family", titleFont.family)
        .attr("font-size", "20px")
        .attr("font-weight", "bold")
        // .attr("text-align", "center")
        .attr("x", (width - margin.left - margin.right) / 2)
        .attr("y", margin.top - margin.bottom - 2)
        .text(title);

      g.append("text")
        .attr("text-anchor", "middle")
        .attr("font-family", titleFont.family)
        .attr("font-size", "20px")
        .attr("x", (width - margin.left - margin.right) / 2)
        .attr("y", margin.top - margin.bottom + 18);
    }

    if (axis_labels) {
      // x-axis labels
      g.append("g")
        .attr("transform", `translate(0, ${plotHeight})`)
        .selectAll("text")
        .data(tickText)
        .join("text")
        .attr("x", (_, i) => i * cellSize + cellSize / 2)
        .attr("y", axisGap)
        .attr("dominant-baseline", "middle")
        .attr("text-anchor", "end")
        .attr("font-family", plotFontMonospace.family)
        .attr("font-size", `${axlabel_fontsize}px`)
        .text((txt) => txt)
        .attr(
          "transform",
          (_, i) =>
            `rotate(${270 + axlabel_xrotation}, ${i * cellSize + cellSize / 2}, ${axisGap})`,
        );

      // y-axis labels
      g.append("g")
        .selectAll("text")
        .data(tickText)
        .join("text")
        .attr("x", -axisGap)
        .attr("y", (_, i) => i * cellSize + cellSize / 2)
        .attr("dominant-baseline", "central")
        .attr("text-anchor", "end")
        .attr("font-family", plotFontMonospace.family)
        .attr("font-size", `${axlabel_fontsize}px`)
        .text((txt) => txt)
        .attr(
          "transform",
          (_, i) =>
            `rotate(${360 + axlabel_yrotation}, ${-axisGap}, ${i * cellSize + cellSize / 2})`,
        );
    }

    if (showscale) {
      const scaleBarX = width - cbarWidth - margin.left - margin.right;

      const gradientGroup = g
        .append("g")
        .attr("transform", `translate(${scaleBarX}, 0)`);
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

      gradientGroup
        .append("rect")
        .attr("width", cbarWidth)
        .attr("height", cbarHeight)
        .attr("fill", `url(#${gradientId})`);

      const axis = d3.axisRight(scale).tickValues(tickValues).tickSize(6);

      g.append("g")
        .attr("transform", `translate(${scaleBarX + cbarWidth}, 0)`)
        .call(axis)
        .call((g) => g.select(".domain").remove())
        .call((g) => g.selectAll(".tick line").attr("stroke-opacity", 0.5))
        .call((g) =>
          g
            .selectAll(".tick text")
            .attr("font-size", "10px")
            .attr("font-family", "Roboto Mono"),
        );
    }
  }, [
    tickValues,
    scale,
    gradientStops,
    cbarHeight,
    cbarWidth,
    svgRef.current,
    data,
    tickText,
    width,
    height,
    cellSpace,
    roundTo,
    showPercentIdentities,
    axlabel_fontsize,
    axlabel_xrotation,
    axlabel_yrotation,
    titleFont,
    showTitles,
    title,
    showscale,
    axis_labels,
    margin,
  ]);

  return (
    <svg
      xmlns="http://www.w3.org/2000/svg"
      xmlnsXlink="http://www.w3.org/1999/xlink"
      style={{ overflow: "visible", background: "#fff" }}
      ref={svgRef}
      width={"100%"}
      height={"100%"}
    />
  );
};
