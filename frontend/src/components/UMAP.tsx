import * as d3 from "d3";
import React, { useEffect, useRef, useState } from "react";
import type { DocState, SetDocState, UpdateDocState } from "../appState";
import type { UMAPData } from "../plotTypes";
import { UMAPSidebar } from "./UMAPSidebar";

interface UMAPProps {
  docState: DocState;
  setDocState: SetDocState;
  updateDocState: UpdateDocState;
  leftSidebarCollapsed: boolean;
}

interface PlotDimensions {
  width: number;
  height: number;
  margin: { top: number; right: number; bottom: number; left: number };
}

interface ExtendedUMAPData extends UMAPData {
  params?: {
    n_neighbors: number;
    min_dist: number;
  };
}

interface SVGElementWithTooltip extends SVGSVGElement {
  __tooltip?: d3.Selection<HTMLDivElement, unknown, HTMLElement, any>;
}

export const UMAP: React.FC<UMAPProps> = ({
  docState,
  setDocState,
  leftSidebarCollapsed,
}) => {
  const containerRef = useRef<HTMLDivElement>(null);
  const svgRef = useRef<SVGElementWithTooltip | null>(null);
  const [umapData, setUmapData] = useState<ExtendedUMAPData | null>(null);
  const [loading, setLoading] = useState(false);
  const [dimensions, setDimensions] = useState({ width: 800, height: 600 });

  // Update settings
  const updateSettings = React.useCallback(
    (values: Partial<DocState["umap"]>) =>
      setDocState((prev) => ({
        ...prev,
        umap: {
          ...prev.umap,
          ...values,
        },
      })),
    [setDocState],
  );

  // Fetch UMAP data when projection parameters change
  useEffect(() => {
    if (!docState.id) return;

    // Only fetch if we don't have data or if parameters changed
    const shouldFetch =
      !umapData ||
      umapData.params?.n_neighbors !== docState.umap.n_neighbors ||
      umapData.params?.min_dist !== docState.umap.min_dist;

    if (!shouldFetch) return;

    setLoading(true);
    window.pywebview.api
      .get_umap_data(docState.id, {
        n_neighbors: docState.umap.n_neighbors,
        min_dist: docState.umap.min_dist,
      })
      .then((response) => {
        setUmapData({
          ...response.data,
          params: {
            n_neighbors: docState.umap.n_neighbors,
            min_dist: docState.umap.min_dist,
          },
        });
      })
      .catch((error) => {
        console.error("Error fetching UMAP data:", error);
      })
      .finally(() => {
        setLoading(false);
      });
  }, [docState.id, docState.umap.n_neighbors, docState.umap.min_dist]);

  // Handle resize
  useEffect(() => {
    const handleResize = () => {
      if (containerRef.current) {
        const { width, height } = containerRef.current.getBoundingClientRect();
        setDimensions({ width, height });
      }
    };

    handleResize();
    window.addEventListener("resize", handleResize);
    return () => window.removeEventListener("resize", handleResize);
  }, []);

  // Calculate plot dimensions
  const getPlotDimensions = (): PlotDimensions => {
    const margin = { top: 40, right: 40, bottom: 50, left: 50 };
    return {
      width: dimensions.width - margin.left - margin.right,
      height: dimensions.height - margin.top - margin.bottom,
      margin,
    };
  };

  // Main render effect
  useEffect(() => {
    if (!umapData || loading || !svgRef.current) return;

    const svg = d3.select(svgRef.current);
    svg.selectAll("*").remove();

    const { width, height, margin } = getPlotDimensions();

    // Main group for the plot
    const mainGroup = svg
      .append("g")
      .attr("transform", `translate(${margin.left},${margin.top})`);

    // Title
    mainGroup
      .append("text")
      .attr("x", width / 2)
      .attr("y", -20)
      .attr("text-anchor", "middle")
      .style("font-size", "18px")
      .style("font-weight", "bold")
      .text("UMAP Projection");

    // Create scales
    const xScale = d3.scaleLinear().domain(umapData.bounds.x).range([0, width]);
    const yScale = d3
      .scaleLinear()
      .domain(umapData.bounds.y)
      .range([height, 0]);

    // Add axes
    mainGroup
      .append("g")
      .attr("transform", `translate(0,${height})`)
      .call(d3.axisBottom(xScale).ticks(5))
      .attr("class", "x-axis");

    mainGroup
      .append("g")
      .call(d3.axisLeft(yScale).ticks(5))
      .attr("class", "y-axis");

    // Create container for plot content
    const plotArea = mainGroup.append("g").attr("class", "plot-area");

    // Create tooltip
    const tooltip = d3
      .select("body")
      .append("div")
      .attr("class", "umap-tooltip")
      .style("position", "absolute")
      .style("visibility", "hidden")
      .style("background", "rgba(0, 0, 0, 0.8)")
      .style("color", "white")
      .style("padding", "8px")
      .style("border-radius", "4px")
      .style("font-size", "12px")
      .style("pointer-events", "none")
      .style("z-index", "1000");

    // Render points
    plotArea
      .selectAll("circle")
      .data(umapData.embedding)
      .enter()
      .append("circle")
      .attr("cx", (d) => xScale(d.x))
      .attr("cy", (d) => yScale(d.y))
      .attr("r", docState.umap.pointSize)
      .attr("fill", "#1f77b4")
      .attr("opacity", docState.umap.opacity)
      .attr("stroke", "white")
      .attr("stroke-width", 0.5)
      .on("mouseover", (_event, d) => {
        tooltip.style("visibility", "visible").html(`ID: ${d.id}`);
      })
      .on("mousemove", (event) => {
        tooltip
          .style("top", `${event.pageY - 10}px`)
          .style("left", `${event.pageX + 10}px`);
      })
      .on("mouseout", () => {
        tooltip.style("visibility", "hidden");
      });

    // Add axis labels
    mainGroup
      .append("text")
      .attr("transform", "rotate(-90)")
      .attr("y", 0 - margin.left)
      .attr("x", 0 - height / 2)
      .attr("dy", "1em")
      .style("text-anchor", "middle")
      .text("UMAP 2");

    mainGroup
      .append("text")
      .attr("transform", `translate(${width / 2}, ${height + margin.bottom})`)
      .style("text-anchor", "middle")
      .text("UMAP 1");

    // Store cleanup function
    const svgNode = svg.node();
    if (svgNode) {
      (svgNode as any).__tooltip = tooltip;
    }

    // Cleanup on unmount
    return () => {
      if (svgRef.current && (svgRef.current as any).__tooltip) {
        (svgRef.current as any).__tooltip.remove();
      }
    };
  }, [umapData, dimensions, docState.umap, loading]);

  return (
    <>
      <div className="app-main umap-container" ref={containerRef}>
        {loading && (
          <div className="app-loader app-main-loader" aria-hidden="true">
            Computing UMAP projection...
          </div>
        )}
        <svg
          ref={svgRef}
          width={dimensions.width}
          height={dimensions.height}
          style={{ maxWidth: "100%", maxHeight: "100%" }}
        />
      </div>
      <UMAPSidebar
        settings={docState.umap}
        updateSettings={updateSettings}
        leftSidebarCollapsed={leftSidebarCollapsed}
      />
    </>
  );
};
