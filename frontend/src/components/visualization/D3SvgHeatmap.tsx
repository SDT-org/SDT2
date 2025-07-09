import * as d3 from "d3";
import React from "react";
import { distinctColor } from "../../colors";
import { plotFontMonospace } from "../../constants";
import { getCellMetrics } from "../../heatmapUtils";
import { useExportSvg } from "../../hooks/useExportSvg";
import type { HeatmapRenderProps } from "../visualization/Heatmap";

export const D3SvgHeatmap = ({
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
  clusterData,
  showLegend,
  showClusterCounts,
  clusterCounts,
}: HeatmapRenderProps) => {
  const svgRef = React.useRef<SVGSVGElement>(null);
  const exportSvg = useExportSvg(
    clusterData ? "clustermap" : "heatmap",
    svgRef,
  );

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

  const rendering = React.useRef(false);
  const onRenderComplete = React.useCallback(() => {
    rendering.current = false;
    const id = setTimeout(() => {
      exportSvg();
    }, 0);
    return () => {
      rendering.current = false;
      clearTimeout(id);
    };
  }, [exportSvg]);

  const renderHeatmap = React.useCallback(
    (ref: React.RefObject<SVGSVGElement>) => {
      if (rendering.current) {
        console.warn("Rendering already in progress");
        return;
      }
      rendering.current = true;

      const start = performance.now();
      const node = ref.current as Element;
      const d3Svg = d3.select(node);
      d3Svg.selectAll("*").remove();
      const plotSize = Math.min(width, height);
      const plotWidth = plotSize - margin.left - margin.right;
      const plotHeight = plotSize - margin.top - margin.bottom;
      const cellSize = plotWidth / tickText.length;
      const cellMetrics = getCellMetrics(cellSize, cellSpace, roundTo + 3);

      const g = d3
        .select(node)
        .append("g")
        .attr("transform", `translate(${margin.left}, ${margin.top})`);

      const axisGap = 5;

      g.selectAll("rect.cell")
        .data(data)
        .join("rect")
        .attr("class", "cell")
        .attr("width", Math.max(cellMetrics.cellSize, 1))
        .attr("height", Math.max(cellMetrics.cellSize, 1))
        .attr("x", (d) => d.x * cellSize + cellMetrics.cellOffset)
        .attr("y", (d) => d.y * cellSize + cellMetrics.cellOffset)
        .attr("fill", (d) => d.backgroundColor);

      if (showPercentIdentities) {
        g.selectAll("text.cell-text")
          .data(data)
          .join("text")
          .attr("class", "cell-text")
          .attr("x", (d) => d.x * cellSize + cellSize / 2)
          .attr("y", (d) => d.y * cellSize + cellSize / 2)
          .attr("dy", ".35em")
          .attr("text-anchor", "middle")
          .attr("font-family", "Roboto Mono")
          .attr("font-size", `${cellMetrics.fontSize}px`)
          .text((d) => d.displayValue)
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
          .on(
            "zoom",
            ({ transform }: d3.D3ZoomEvent<SVGSVGElement, unknown>) => {
              g.attr("transform", transform.toString());
            },
          ),
      );

      if (showTitles) {
        g.append("text")
          .attr("text-anchor", "middle")
          .attr("font-family", titleFont.family)
          .attr("font-size", "20px")
          .attr("font-weight", "bold")
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
        const scaleBarX = Math.min(
          plotSize - cbarWidth - margin.left,
          width - cbarWidth - margin.left - margin.right,
        );

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

      if (clusterData && showLegend) {
        const legendWidth = showClusterCounts ? 120 : 80;
        const cellSize = 10;
        const lineGap = 20;
        const labelGap = 5;
        const columnGap = showClusterCounts ? 30 : 20;
        const positionX = Math.min(
          plotSize - margin.left,
          width - legendWidth * 2 - columnGap - margin.right,
        );

        const uniqueClusters = [...new Set(clusterData.map((i) => i.cluster))]
          .sort((a, b) => a - b)
          .slice(0, 50);

        const legends = g
          .append("g")
          .attr("transform", () => `translate(-${margin.right}, 0)`)
          .selectAll("g")
          .data(uniqueClusters)
          .join("g")
          .attr("transform", (d) => {
            const index = d - 1; // data from d3 is 1-indexed
            const column = index % 2;
            const row = Math.floor(index / 2);
            const itemX = positionX + column * (legendWidth + columnGap);
            const itemY = lineGap * row;
            return `translate(${itemX}, ${itemY})`;
          });

        legends
          .append("rect")
          .attr("width", cellSize)
          .attr("height", cellSize)
          .attr("x", 0)
          .attr("y", 0)
          .attr("fill", (d) => distinctColor(d));

        if (showClusterCounts && clusterCounts) {
          // Cluster label (left-aligned)
          legends
            .append("text")
            .attr("x", cellSize + labelGap)
            .attr("y", cellSize / 2)
            .attr("dy", ".35em")
            .attr("font-family", "Roboto Mono")
            .attr("font-size", "10px")
            .attr("text-anchor", "start")
            .text((d) => {
              // Format cluster number with leading space for single digits
              const clusterStr = d < 10 ? ` ${d}` : `${d}`;
              return `Cluster ${clusterStr}`;
            })
            .attr("fill", "black");

          // Count (right-aligned)
          legends
            .append("text")
            .attr("x", legendWidth)
            .attr("y", cellSize / 2)
            .attr("dy", ".35em")
            .attr("font-family", "Roboto Mono")
            .attr("font-size", "10px")
            .attr("text-anchor", "end")
            .text((d) => {
              if (clusterCounts[d]) {
                // Add leading space before bracket for single digits
                const countDisplay =
                  clusterCounts[d] < 10
                    ? ` [${clusterCounts[d]}]`
                    : `[${clusterCounts[d]}]`;
                return countDisplay;
              }
              return "";
            })
            .attr("fill", "black");
        } else {
          legends
            .append("text")
            .attr("x", cellSize + labelGap)
            .attr("y", cellSize / 2)
            .attr("dy", ".35em")
            .attr("font-family", "Roboto Mono")
            .attr("font-size", "10px")
            .text((d) => {
              // Format cluster number with leading space for single digits
              const clusterStr = d < 10 ? ` ${d}` : `${d}`;
              return `Cluster ${clusterStr}`;
            })
            .attr("fill", "black");
        }
      }

      const end = performance.now();
      const duration = end - start;
      console.log(`Heatmap rendering took ${duration.toFixed(2)} milliseconds`);

      onRenderComplete();
    },
    [
      onRenderComplete,
      tickValues,
      scale,
      gradientStops,
      cbarHeight,
      cbarWidth,
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
      clusterData,
      showLegend,
      showClusterCounts,
      clusterCounts,
    ],
  );

  React.useEffect(() => {
    if (!(svgRef.current && height && width)) return;
    renderHeatmap(svgRef);
  }, [renderHeatmap, height, width]);

  return (
    <svg
      xmlns="http://www.w3.org/2000/svg"
      xmlnsXlink="http://www.w3.org/1999/xlink"
      id={"heatmap-svg"}
      style={{ overflow: "visible", background: "#fff" }}
      ref={svgRef}
      width={"100%"}
      height={"100%"}
    />
  );
};
