import * as d3 from "d3";
import type React from "react";
import { useEffect, useMemo, useRef } from "react";
import type { ColorScaleArray } from "../colorScales";

//entirely based on https://observablehq.com/@d3/color-legend
interface ColorLegendProps {
  colorScale: ColorScaleArray;
  minVal: number;
  maxVal: number;
  showscale?: boolean;
  position?: { x: number; y: number };
  title?: string;
  ticks?: number;
  cbarHeight?: number;
  cbarWidth?: number;
}

export const ColorLegend: React.FC<ColorLegendProps> = ({
  colorScale,
  minVal,
  maxVal,
  showscale = true,
  cbarWidth = 60,
  cbarHeight = 200,
  position = { x: 0, y: 0 },
  title = "",
  ticks = 5,
}) => {
  const svgRef = useRef<SVGSVGElement>(null);

  const dimensions = useMemo(
    () => ({ cbarWidth, cbarHeight }),
    [cbarWidth, cbarHeight],
  );

  useEffect(() => {
    if (!svgRef.current || !showscale) return;

    const svg = d3.select(svgRef.current);
    svg.selectAll("*").remove();

    const margin = { top: 20, right: 30, bottom: 20, left: 0 };
    const g = svg
      .append("g")
      .attr("transform", `translate(${margin.left},${margin.top})`);

    const defs = svg.append("defs");
    const gradientId = `legend-gradient-${Math.random().toString(36).substr(2, 9)}`;
    const gradient = defs
      .append("linearGradient")
      .attr("id", gradientId)
      .attr("x1", "0%")
      .attr("y1", "100%")
      .attr("x2", "0%")
      .attr("y2", "0%");

    for (const [stop, color] of colorScale) {
      gradient
        .append("stop")
        .attr("offset", `${stop * 100}%`)
        .attr("stop-color", color);
    }
    
    const legendHeight = dimensions.cbarHeight - margin.top - margin.bottom;
    const legendWidth = 50; //padding

    g.append("rect")
      .attr("width", legendWidth)
      .attr("height", legendHeight)
      .style("fill", `url(#${gradientId})`);

    const scale = d3
      .scaleLinear()
      .domain([maxVal, minVal]) // Reversed domain for vertical axis ticks
      .range([0, legendHeight]);

    const axis = d3.axisRight(scale).ticks(ticks).tickSize(6);

    g.append("g")
      .attr("transform", `translate(${legendWidth},0)`)
      .call(axis)
      .call((g) => g.select(".domain").remove())
      .call((g) => g.selectAll(".tick line").attr("stroke-opacity", 0.5))
      .call((g) => g.selectAll(".tick text").attr("font-size", "10px"));

    if (title) {
      g.append("text")
        .attr("class", "legend-title")
        .attr("x", -legendHeight / 2)
        .attr("y", -margin.left)
        .attr("transform", "rotate(-90)")
        .attr("text-anchor", "middle")
        .attr("font-size", "12px")
        .attr("font-weight", "bold")
        .text(title);
    }
  }, [colorScale, minVal, maxVal, showscale, title, ticks, dimensions]);

  if (!showscale) return null;

  return (
    <div
      style={{
        position: "absolute",
        top: position.y,
        right: position.x,
        background: "transparent",
        width: `${dimensions.cbarWidth}px`,
        height: `${dimensions.cbarHeight}px`,
      }}
    >
      <svg
        ref={svgRef}
        width={dimensions.cbarWidth}
        height={dimensions.cbarHeight}
        style={{ overflow: "visible" }}
      />
    </div>
  );
};
