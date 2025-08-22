import * as d3 from "d3";
import type React from "react";
import { useEffect } from "react";
import { distinctColor } from "../../../colors";
import { getMetadataTooltip } from "../../../hooks/useMetadataColors";
import type { UMAPPoint } from "../../../plotTypes";

interface D3SvgUMAPProps {
  umapData: {
    embedding: UMAPPoint[];
    bounds: { x: [number, number]; y: [number, number] };
    clusterStats?: {
      total_clusters: number;
      noise_points: number;
      largest_cluster_size: number;
      smallest_cluster_size: number;
    };
  };
  dimensions: { width: number; height: number };
  pointSize: number;
  opacity: number;
  colorBy: string;
  colorByCluster: boolean;
  selectedMetadataColumn?: string;
  metadataValues?: Record<string, string | number | boolean | null>;
  selectedPoints: string[];
  selectionActive: boolean;
  polygon: [number, number][];
  svgRef: React.RefObject<SVGSVGElement>;
  getMetadataColor: (id: string) => string;
  setPolygon: (polygon: [number, number][]) => void;
  setSelectedPoints: (points: string[]) => void;
  setSelectionActive: (active: boolean) => void;
  setSelectionSummary: (
    summary: {
      summary:
        | Record<string, number>
        | {
            mean: number;
            median: number;
            min: number;
            max: number;
            [key: string]: number;
          };
      column_type: "numeric" | "categorical";
      total_selected: number;
    } | null,
  ) => void;
  docId: string;
  error: string | null;
  loading: boolean;
}

export const D3SvgUMAP: React.FC<D3SvgUMAPProps> = ({
  umapData,
  dimensions,
  pointSize,
  opacity,
  colorBy,
  colorByCluster,
  selectedMetadataColumn,
  metadataValues,
  selectedPoints,
  selectionActive,
  polygon,
  svgRef,
  getMetadataColor,
  setPolygon,
  setSelectedPoints,
  setSelectionActive,
  setSelectionSummary,
  docId,
  error,
  loading,
}) => {
  // Main render effect - only renders the base visualization
  // biome-ignore lint/correctness/useExhaustiveDependencies: Ref is stable and D3 manages its own DOM
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

    // Handle brush selection for metadata summary
    const handleBrush = (event: d3.D3BrushEvent<unknown>) => {
      const { selection } = event;
      const allPoints = mainGroup.selectAll(".umap-point");

      if (!selection) {
        if (!selectionActive) {
          allPoints.classed("selected", false);
        }
        return;
      }

      const [[x0, y0], [x1, y1]] = selection as [
        [number, number],
        [number, number],
      ];

      // Cast the datum to UMAPPoint but use proper typing for D3
      allPoints.classed("selected", (d) => {
        const point = d as unknown as UMAPPoint;
        const px = currentXScale(point.x);
        const py = currentYScale(point.y);
        return x0 <= px && px <= x1 && y0 <= py && py <= y1;
      });
    };

    const handleBrushEnd = async (event: d3.D3BrushEvent<unknown>) => {
      const { selection } = event;
      if (!selection) {
        // Only clear if user clicked outside without dragging
        if (!selectionActive) {
          setSelectionSummary(null);
          setSelectedPoints([]);
          setSelectionActive(false);
          mainGroup.selectAll(".umap-point").classed("selected", false);
        }
        return;
      }

      const [[x0, y0], [x1, y1]] = selection as [
        [number, number],
        [number, number],
      ];

      const pointsInSelection = umapData.embedding.filter((d) => {
        const px = currentXScale(d.x);
        const py = currentYScale(d.y);
        return x0 <= px && px <= x1 && y0 <= py && py <= y1;
      });

      const selectedIds = pointsInSelection.map((d) => d.id);

      if (selectedIds.length > 0) {
        setSelectedPoints(selectedIds);
        setSelectionActive(true);

        if (selectedMetadataColumn) {
          try {
            const summary =
              await window.pywebview.api.data.get_metadata_summary_for_selection(
                docId,
                selectedIds,
                selectedMetadataColumn,
              );
            setSelectionSummary(summary);
          } catch (err) {
            console.error("Error fetching selection summary:", err);
          }
        }
      }

      // Clear the brush selection but maintain the visual highlight
      brush.move(mainGroup.select(".brush"), null);
    };

    // Render points with initial positions
    plotArea
      .selectAll("circle")
      .data(umapData.embedding)
      .enter()
      .append("circle")
      .attr("cx", (d) => xScale(d.x))
      .attr("cy", (d) => yScale(d.y))
      .attr("r", pointSize)
      .attr("fill", (d) => {
        if (colorBy === "metadata") {
          return getMetadataColor(d.id);
        }

        // Default cluster coloring
        if (!colorByCluster || !d.cluster) {
          return "#1f77b4"; // Default blue
        }
        const clusterId = d.cluster ?? 0;
        return clusterId === 0 ? "#cccccc" : distinctColor(clusterId); // Gray for noise points
      })
      .attr("opacity", opacity)
      .attr("stroke", "white")
      .attr("stroke-width", 0.5)
      .attr("class", "umap-point") // Add class for easy selection
      .on("mouseover", (_event, d) => {
        tooltip.style("visibility", "visible");
        let content = `<strong>ID:</strong> ${d.id}`;
        if (d.cluster) {
          content += `<br><strong>Cluster:</strong> ${d.cluster}`;
        }
        if (colorBy === "metadata" && metadataValues) {
          const metadataValue = metadataValues[d.id];
          if (
            metadataValue !== null &&
            metadataValue !== undefined &&
            typeof metadataValue !== "boolean"
          ) {
            content += getMetadataTooltip(
              d.id,
              selectedMetadataColumn,
              metadataValue as string | number,
            );
          }
        }
        tooltip.html(content);
      })
      .on("mousemove", (event) => {
        tooltip
          .style("top", `${event.pageY - 10}px`)
          .style("left", `${event.pageX + 10}px`);
      })
      .on("mouseout", () => {
        tooltip.style("visibility", "hidden");
      });

    // Add brush - properly initialized with explicit handlers
    const brush = d3
      .brush()
      .extent([
        [0, 0],
        [width, height],
      ])
      .filter((event) => event.button === 2) // Only brush with right-click
      .on("start", handleBrush)
      .on("brush", handleBrush)
      .on("end", handleBrushEnd);

    mainGroup
      .append("g")
      .attr("class", "brush")
      .attr("aria-label", "Selection brush area")
      .call(brush);

    // Prevent context menu on right-click to allow for other interactions
    svg.on("contextmenu", (event) => event.preventDefault());

    // Add polygon selection with Ctrl+right-click
    svg.on("mousedown", (event) => {
      // Only add points with Ctrl+right-click
      if (event.ctrlKey && event.button === 2) {
        // Prevent default to avoid triggering other handlers
        event.preventDefault();

        const [mx, my] = d3.pointer(event);
        const newPolygon = [
          ...polygon,
          [mx - margin.left, my - margin.top],
        ] as [number, number][];
        setPolygon(newPolygon);
      }
    });

    // Add help text for selection methods
    mainGroup
      .append("text")
      .attr("x", 10)
      .attr("y", -10)
      .style("font-size", "12px")
      .style("fill", "#666")
      .text(
        "Right-click drag to brush select. Ctrl+Right-click to draw polygon.",
      );

    // Show polygon path if points exist
    if (polygon.length > 0) {
      mainGroup
        .append("path")
        .datum(polygon)
        .attr("class", "selection-polygon")
        .attr("d", d3.line())
        .attr("stroke", "black")
        .attr("stroke-width", 2)
        .attr("fill", "rgba(0,0,0,0.1)");
    }

    // Show red start point if polygon has at least 3 points
    if (polygon.length > 2) {
      const firstPoint = mainGroup
        .selectAll<SVGCircleElement, [number, number]>(".polygon-start-point")
        .data([polygon[0]]);

      firstPoint
        .enter()
        .append("circle")
        .attr("class", "polygon-start-point")
        .attr("r", 5)
        .attr("fill", "red")
        .merge(firstPoint)
        .attr("cx", (d) => (d ? d[0] : null))
        .attr("cy", (d) => (d ? d[1] : null))
        .on("mousedown", async (event) => {
          // Only close polygon with Ctrl+right-click (consistent with adding points)
          if (event.ctrlKey && event.button === 2) {
            // Prevent default to avoid unexpected behavior
            event.preventDefault();

            const pointsInPolygon = umapData.embedding.filter((d) => {
              const point = [currentXScale(d.x), currentYScale(d.y)] as [
                number,
                number,
              ];
              return d3.polygonContains(polygon, point);
            });
            const selectedIds = pointsInPolygon.map((d) => d.id);

            if (selectedIds.length > 0) {
              setSelectedPoints(selectedIds);
              setSelectionActive(true);

              // For metadata if available - always fetch summary regardless of current coloring mode
              if (selectedMetadataColumn) {
                try {
                  const summary =
                    await window.pywebview.api.data.get_metadata_summary_for_selection(
                      docId,
                      selectedIds,
                      selectedMetadataColumn,
                    );
                  setSelectionSummary(summary);
                } catch (err) {
                  console.error("Error fetching selection summary:", err);
                }
              }
            }

            setPolygon([]);
          }
        });
    }

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
        xAxisGroup.call(d3.axisBottom(currentXScale).ticks(5));
        yAxisGroup.call(d3.axisLeft(currentYScale).ticks(5));

        // Update points positions
        svg
          .selectAll<SVGCircleElement, UMAPPoint>(".umap-point")
          .attr("cx", (d) => currentXScale(d.x))
          .attr("cy", (d) => currentYScale(d.y));
      });

    // Apply zoom to SVG
    svg.call(zoom as unknown as d3.ZoomBehavior<SVGSVGElement, unknown>);

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
        // d3 typing is complex between transitions and selections
        svg
          .transition()
          .duration(750)
          // @ts-ignore d3 typings don't properly handle zoom transforms on transitions
          .call(zoom.transform, d3.zoomIdentity);
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

    // Restore selection state if there are selected points
    if (selectionActive && selectedPoints.length > 0) {
      plotArea
        .selectAll<SVGCircleElement, UMAPPoint>(".umap-point")
        .classed("selected", (d) => selectedPoints.includes(d.id));
    }

    // Cleanup on unmount
    return () => {
      tooltip.remove();
    };
  }, [
    umapData,
    dimensions,
    loading,
    error,
    docId,
    colorBy,
    colorByCluster,
    selectedMetadataColumn,
    metadataValues,
    pointSize,
    opacity,
    getMetadataColor,
    selectedPoints,
    selectionActive,
    polygon,
    setPolygon,
    setSelectedPoints,
    setSelectionActive,
    setSelectionSummary,
  ]);

  // Update cluster stats panel separately
  // biome-ignore lint/correctness/useExhaustiveDependencies: Ref is stable and D3 manages its own DOM
  useEffect(() => {
    if (!svgRef.current || !umapData) return;

    const svg = d3.select(svgRef.current);
    const mainGroup = svg.select("g");

    // Remove existing stats group
    mainGroup.select(".cluster-stats").remove();

    const clusterStats = umapData.clusterStats;
    const width = dimensions.width - 90; // margin.left + margin.right

    // Add cluster stats panel
    if (clusterStats && colorByCluster) {
      const statsGroup = mainGroup
        .append("g")
        .attr("class", "cluster-stats")
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
  }, [umapData, dimensions, colorByCluster]);

  // Return an empty div if we're not rendering
  if (loading || error || !umapData) {
    return null;
  }

  return (
    <svg
      ref={svgRef}
      width={dimensions.width}
      height={dimensions.height}
      style={{ display: loading ? "none" : "block" }}
    />
  );
};
