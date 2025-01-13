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
  tempHeatmapComponent: "canvas" | "svg" | "plotly"; 
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
  tempHeatmapComponent = "svg", 
  cbarWidth = 60,
  cbarHeight = 200,
  position = { x: 0, y: 0 },
  title = "",
  ticks = 5,
}) => {
  const svgRef = useRef<SVGSVGElement>(null);
  const canvasRef = useRef<HTMLCanvasElement>(null);

  const dimensions = useMemo(
    () => ({ cbarWidth, cbarHeight }),
    [cbarWidth, cbarHeight]
  );

  // Shared logic to compute gradient stops, scale, and ticks
  const gradientStops = useMemo(() => colorScale.map(([stop, color]) => ({ stop, color })), [colorScale]);
  const scale = useMemo(
    () => d3.scaleLinear().domain([maxVal, minVal]).range([0, dimensions.cbarHeight]),
    [minVal, maxVal, dimensions.cbarHeight]
  );
  const tickValues = useMemo(() => scale.ticks(ticks), [scale, ticks]);

  useEffect(() => {
    if (!showscale) return;

    if (tempHeatmapComponent === "svg" && svgRef.current) {
      const svg = d3.select(svgRef.current);
      svg.selectAll("*").remove();

      const margin = { top: 20, right: 30, bottom: 20, left: 0 };
      const g = svg
        .append("g")
        .attr("transform", `translate(${margin.left},${margin.top})`);

      //  gradient
      const defs = svg.append("defs");
      const gradientId = `legend-gradient-${Math.random().toString(36).substr(2, 9)}`;
      const gradient = defs
        .append("linearGradient")
        .attr("id", gradientId)
        .attr("x1", "0%")
        .attr("y1", "0%")
        .attr("x2", "0%")
        .attr("y2", "100%");

      gradientStops.forEach(({ stop, color }) => {
        gradient
          .append("stop")
          .attr("offset", `${stop * 100}%`)
          .attr("stop-color", color);
      });

      // Add the rectangle
      const legendHeight = dimensions.cbarHeight - margin.top - margin.bottom;
      const legendWidth = dimensions.cbarWidth;

      g.append("rect")
        .attr("width", legendWidth)
        .attr("height", legendHeight)
        .attr("fill", `url(#${gradientId})`);

      // Add the axis
      const axis = d3.axisRight(scale).tickValues(tickValues).tickSize(6);
      g.append("g")
        .attr("transform", `translate(${legendWidth},0)`)
        .call(axis)
        .call((g) => g.select(".domain").remove())
        .call((g) => g.selectAll(".tick line").attr("stroke-opacity", 0.5))
        .call((g) => g.selectAll(".tick text").attr("font-size", "10px"));

      // Add title
      if (title) {
        g.append("text")
          .attr("class", "legend-title")
          .attr("x", -legendHeight / 2)
          .attr("y", -margin.left - 10)
          .attr("transform", "rotate(-90)")
          .attr("text-anchor", "middle")
          .attr("font-size", "12px")
          .attr("font-weight", "bold")
          .text(title);
      }
    } else if (tempHeatmapComponent === "canvas" && canvasRef.current) {
      const canvas = canvasRef.current;
      const context = canvas.getContext("2d");
      if (!context) return;

      // Clear canvas
      context.clearRect(0, 0, canvas.width, canvas.height);

      // Create the gradient
      const gradient = context.createLinearGradient(0, 0, 0, cbarHeight);
      gradientStops.forEach(({ stop, color }) => {
        gradient.addColorStop(stop, color);
      });

      // Fill the gradient
      context.fillStyle = gradient;
      context.fillRect(0, 0, cbarWidth, cbarHeight);

      // Draw ticks
      context.fillStyle = "#000";
      context.font = "10px sans-serif";
      context.textAlign = "left";

      tickValues.forEach((tick) => {
        const y = scale(tick);
        context.fillText(tick.toFixed(2), cbarWidth + 5, y + 3);
        context.beginPath();
        context.moveTo(cbarWidth, y);
        context.lineTo(cbarWidth - 6, y);
        context.strokeStyle = "#000";
        context.stroke();
      });
    }
  }, [gradientStops, scale, tickValues, tempHeatmapComponent, showscale, dimensions, title]);

  if (!showscale) return null;

  return (
    <div
      style={{
        position: "absolute",
        top: position.y,
        left: position.x,
        width: `${dimensions.cbarWidth}px`,
        height: `${dimensions.cbarHeight}px`,
      }}
    >
      {tempHeatmapComponent === "svg" ? (
        <svg
          ref={svgRef}
          width={dimensions.cbarWidth}
          height={dimensions.cbarHeight}
          style={{ overflow: "visible" }}
        />
      ) : (
        <canvas
          ref={canvasRef}
          width={dimensions.cbarWidth}
          height={dimensions.cbarHeight}
          style={{ overflow: "visible" }}
        />
      )}
    </div>
  );
};
