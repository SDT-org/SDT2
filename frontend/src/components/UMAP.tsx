import * as d3 from "d3";
import React, { useEffect, useRef, useState } from "react";
import type { DocState, SetDocState, UpdateDocState } from "../appState";
import { distinctColor } from "../colors";
import type { UMAPData, UMAPPoint } from "../plotTypes";
import { UMAPSidebar } from "./UMAPSidebar";

interface UMAPProps {
  docState: DocState;
  setDocState: SetDocState;
  updateDocState: UpdateDocState;
  leftSidebarCollapsed: boolean;
}

interface ExtendedUMAPData extends UMAPData {
  params?: {
    n_neighbors: number;
    min_dist: number;
  };
  clusterStats?: {
    total_clusters: number;
    noise_points: number;
    largest_cluster_size: number;
    smallest_cluster_size: number;
  };
}

export const UMAP: React.FC<UMAPProps> = ({
  docState,
  setDocState,
  leftSidebarCollapsed,
}) => {
  const containerRef = useRef<HTMLDivElement>(null);
  const svgRef = useRef<SVGSVGElement | null>(null);
  const [umapData, setUmapData] = useState<ExtendedUMAPData | null>(null);
  const [loading, setLoading] = useState(false);
  const [dimensions, setDimensions] = useState({ width: 800, height: 600 });
  const [error, setError] = useState<string | null>(null);

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

  // Fetch UMAP data when parameters change
  useEffect(() => {
    if (!docState.id) return;

    const fetchData = async () => {
      setLoading(true);
      setError(null);

      try {
        const response = await window.pywebview.api.data.get_umap_data(
          docState.id,
          {
            n_neighbors: docState.umap.n_neighbors,
            min_dist: docState.umap.min_dist,
            threshold: docState.umap.minClusterSize, // Backend interprets this as min_cluster_size
            methods: [`hdbscan-${docState.umap.clusterEpsilon}`], // Backend extracts epsilon from method string
          },
        );

        // Debug logging
        console.log("UMAP data received:", response.data);
        if (response.data.embedding && response.data.embedding.length > 0) {
          console.log(
            "First 5 points cluster data:",
            response.data.embedding.slice(0, 5),
          );
        }

        setUmapData({
          ...response.data,
          params: {
            n_neighbors: docState.umap.n_neighbors,
            min_dist: docState.umap.min_dist,
          },
        });
      } catch (err) {
        console.error("Error fetching UMAP data:", err);
        const errorMessage =
          err instanceof Error ? err.message : "Failed to load UMAP data";
        setError(errorMessage);
      } finally {
        setLoading(false);
      }
    };

    fetchData();
  }, [
    docState.id,
    docState.umap.n_neighbors,
    docState.umap.min_dist,
    docState.umap.minClusterSize,
    docState.umap.clusterEpsilon,
  ]);

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

  // Main render effect
  useEffect(() => {
    if (!umapData || loading || !svgRef.current || error) return;

    const svg = d3.select(svgRef.current);
    svg.selectAll("*").remove();

    // Calculate plot dimensions
    const margin = { top: 60, right: 40, bottom: 50, left: 50 };
    const width = dimensions.width - margin.left - margin.right;
    const height = dimensions.height - margin.top - margin.bottom;

    // Main group for the plot
    const mainGroup = svg
      .append("g")
      .attr("transform", `translate(${margin.left},${margin.top})`);

    // Title with cluster stats
    const clusterStats = umapData.clusterStats;
    const title = clusterStats
      ? `UMAP Projection - HDBSCAN (${clusterStats.total_clusters} clusters, ${clusterStats.noise_points} noise points)`
      : "UMAP Projection";

    mainGroup
      .append("text")
      .attr("x", width / 2)
      .attr("y", -30)
      .attr("text-anchor", "middle")
      .style("font-size", "18px")
      .style("font-weight", "bold")
      .text(title);

    // Create scales with copies for zoom
    const xScale = d3.scaleLinear().domain(umapData.bounds.x).range([0, width]);
    const yScale = d3
      .scaleLinear()
      .domain(umapData.bounds.y)
      .range([height, 0]);

    let currentXScale = xScale.copy();
    let currentYScale = yScale.copy();

    // Add axes
    const xAxisGroup = mainGroup
      .append("g")
      .attr("transform", `translate(0,${height})`)
      .attr("class", "x-axis");

    const yAxisGroup = mainGroup.append("g").attr("class", "y-axis");

    const xAxis = d3.axisBottom(currentXScale).ticks(5);
    const yAxis = d3.axisLeft(currentYScale).ticks(5);

    xAxisGroup.call(xAxis);
    yAxisGroup.call(yAxis);

    // Add clip path to prevent drawing outside plot area
    svg
      .append("defs")
      .append("clipPath")
      .attr("id", "plot-clip")
      .append("rect")
      .attr("x", margin.left)
      .attr("y", margin.top)
      .attr("width", width)
      .attr("height", height);

    // Create container for plot content
    const plotArea = mainGroup
      .append("g")
      .attr("class", "plot-area")
      .attr("clip-path", "url(#plot-clip)");

    // Create tooltip
    const tooltip = d3
      .select("body")
      .append("div")
      .attr("class", "umap-tooltip")
      .style("position", "absolute")
      .style("visibility", "hidden")
      .style("background", "rgba(0, 0, 0, 0.9)")
      .style("color", "white")
      .style("padding", "10px")
      .style("border-radius", "4px")
      .style("font-size", "12px")
      .style("pointer-events", "none")
      .style("z-index", "1000");

    // Get color function based on settings
    const getPointColor = (d: UMAPPoint): string => {
      if (!docState.umap.colorByCluster || !d.cluster) {
        return "#1f77b4"; // Default blue
      }
      const clusterId = d.cluster ?? 0;
      return clusterId === 0 ? "#cccccc" : distinctColor(clusterId); // Gray for noise points
    };

    // Render points with initial positions
    const points = plotArea
      .selectAll("circle")
      .data(umapData.embedding)
      .enter()
      .append("circle")
      .attr("cx", (d) => xScale(d.x))
      .attr("cy", (d) => yScale(d.y))
      .attr("r", docState.umap.pointSize)
      .attr("fill", getPointColor)
      .attr("opacity", docState.umap.opacity)
      .attr("stroke", "white")
      .attr("stroke-width", 0.5)
      .on("mouseover", (_event, d) => {
        let tooltipContent = `<strong>ID:</strong> ${d.id}`;
        if (d.cluster !== undefined) {
          const clusterLabel = d.cluster === 0 ? "Noise" : d.cluster;
          tooltipContent += `<br><strong>Cluster:</strong> ${clusterLabel}`;
        }
        tooltip.style("visibility", "visible").html(tooltipContent);
      })
      .on("mousemove", (event) => {
        tooltip
          .style("top", `${event.pageY - 10}px`)
          .style("left", `${event.pageX + 10}px`);
      })
      .on("mouseout", () => {
        tooltip.style("visibility", "hidden");
      });

    // Add zoom behavior
    const zoom = d3
      .zoom<SVGSVGElement, unknown>()
      .scaleExtent([0.5, 10])
      .extent([
        [0, 0],
        [dimensions.width, dimensions.height],
      ])
      .on("zoom", (event) => {
        const { transform } = event;

        // Update scales
        currentXScale = transform.rescaleX(xScale);
        currentYScale = transform.rescaleY(yScale);

        // Update axes
        xAxisGroup.call(xAxis.scale(currentXScale));
        yAxisGroup.call(yAxis.scale(currentYScale));

        // Update points positions
        points
          .attr("cx", (d) => currentXScale(d.x))
          .attr("cy", (d) => currentYScale(d.y));
      });

    // Apply zoom to SVG
    svg.call(zoom);

    // Add reset button
    mainGroup
      .append("text")
      .attr("x", width - 60)
      .attr("y", -30)
      .style("cursor", "pointer")
      .style("font-size", "14px")
      .style("text-decoration", "underline")
      .style("fill", "#1f77b4")
      .text("Reset Zoom")
      .on("click", () => {
        svg.transition().duration(750).call(zoom.transform, d3.zoomIdentity);
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

    // Add cluster stats panel
    if (clusterStats && docState.umap.colorByCluster) {
      const statsGroup = mainGroup
        .append("g")
        .attr("transform", `translate(${width - 200}, 0)`);

      statsGroup
        .append("text")
        .attr("x", 0)
        .attr("y", 0)
        .style("font-size", "14px")
        .style("font-weight", "bold")
        .text("Cluster Statistics");

      const stats = [
        `Total clusters: ${clusterStats.total_clusters}`,
        `Noise points: ${clusterStats.noise_points}`,
        `Largest cluster: ${clusterStats.largest_cluster_size}`,
        `Smallest cluster: ${clusterStats.smallest_cluster_size}`,
      ];

      stats.forEach((stat, i) => {
        statsGroup
          .append("text")
          .attr("x", 0)
          .attr("y", 20 + i * 18)
          .style("font-size", "12px")
          .text(stat);
      });
    }

    // Cleanup on unmount
    return () => {
      tooltip.remove();
    };
  }, [umapData, dimensions, docState.umap, loading, error]);

  return (
    <>
      <div className="app-main umap-container" ref={containerRef}>
        {loading && (
          <>
            <div className="app-loader app-main-loader" aria-hidden="true" />
            <div className="app-loader-status">
              Computing UMAP projection...
            </div>
          </>
        )}
        {error && (
          <div
            className="error-message"
            style={{ padding: "20px", color: "red" }}
          >
            {error}
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
