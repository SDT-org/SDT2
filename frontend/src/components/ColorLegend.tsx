import * as d3 from "d3";
import type React from "react";
import { useEffect, useMemo, useRef } from "react";
import type { ColorScaleArray } from "../colorScales";

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
  cbarWidth = 0,
  cbarHeight = 200,
  position = { x: 0, y: 0 },
  title = "",
  ticks = 5,
}) => {
  const svgRef = useRef<SVGSVGElement>(null);
  const canvasRef = useRef<HTMLCanvasElement>(null);
  const dimensions = useMemo(
    () => ({ cbarWidth, cbarHeight }),
    [cbarWidth, cbarHeight],
  );
  const gradientStops = useMemo(
    () => colorScale.map(([stop, color]) => ({ stop, color })),
    [colorScale],
  );
  const scale = useMemo(
    () =>
      d3
        .scaleLinear()
        .domain([maxVal, minVal])
        .range([0, dimensions.cbarHeight]),
    [minVal, maxVal, dimensions.cbarHeight],
  );
  const tickValues = useMemo(() => scale.ticks(ticks), [scale, ticks]);

  useEffect(() => {
    if (!showscale) return;

    if (tempHeatmapComponent === "svg" && svgRef.current) {
      const svg = d3.select(svgRef.current);
      svg.selectAll("*").remove();

      const margin = { top: 20, right: 30, bottom: 20, left: 10 };
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

      for (const { stop, color } of gradientStops) {
        gradient
          .append("stop")
          .attr("offset", `${stop * 100}%`)
          .attr("stop-color", color);
      }

      const legendHeight = dimensions.cbarHeight - margin.top - margin.bottom;
      const legendWidth = dimensions.cbarWidth - margin.left - margin.right;

      g.append("rect")
        .attr("width", legendWidth)
        .attr("height", legendHeight)
        .attr("fill", `url(#${gradientId})`);

      const axis = d3.axisRight(scale).tickValues(tickValues).tickSize(6);

      g.append("g")
        .attr("transform", `translate(${legendWidth},0)`)
        .call(axis)
        .call((g) => g.select(".domain").remove())
        .call((g) => g.selectAll(".tick line").attr("stroke-opacity", 0.5))
        .call((g) =>
          g
            .selectAll(".tick text")
            .attr("font-size", "10px")
            .attr("dx", "0.5em"),
        );

      if (title) {
        g.append("text")
          .attr("x", -legendHeight / 2)
          .attr("y", -margin.left - 5)
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

      context.clearRect(0, 0, canvas.width, canvas.height);

      const marginLeft = 10;
      const marginRight = 40;
      const marginTop = 10;
      const barHeight = cbarHeight - marginTop * 2;

      const gradient = context.createLinearGradient(
        0,
        barHeight + marginTop,
        0,
        marginTop,
      );
      for (const { stop, color } of gradientStops) {
        gradient.addColorStop(stop, color);
      }

      const barWidth = cbarWidth - marginLeft - marginRight;
      context.fillStyle = gradient;
      context.fillRect(marginLeft, marginTop, barWidth, barHeight);

      const adjustedScale = scale
        .copy()
        .range([marginTop, barHeight + marginTop]);

      context.fillStyle = "#000";
      context.font = "10px 'Roboto Mono'";
      context.textAlign = "left";
      context.textBaseline = "middle";

      for (const tick of tickValues) {
        const y = adjustedScale(tick);
        context.fillText(tick.toFixed(2), marginLeft + barWidth + 5, y);
        context.beginPath();
        context.moveTo(marginLeft + barWidth, y);
        context.lineTo(marginLeft + barWidth + 6, y);
        context.strokeStyle = "#000";
        context.stroke();
      }
    }
  }, [
    gradientStops,
    scale,
    tickValues,
    tempHeatmapComponent,
    showscale,
    dimensions,
    title,
    cbarWidth,
    cbarHeight,
  ]);

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
          style={{ overflow: "hidden" }}
        />
      )}
    </div>
  );
};
