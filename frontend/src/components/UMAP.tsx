import * as d3 from "d3";
import React, { useEffect, useRef, useState } from "react";
import type { DocState, SetDocState, UpdateDocState } from "../appState";
import { distinctColor } from "../colors";
import { useExportSvg } from "../hooks/useExportSvg";
import {
  getMetadataTooltip,
  useMetadataColors,
} from "../hooks/useMetadataColors";
import type { UMAPPoint } from "../plotTypes";
import { UMAPSidebar } from "./UMAPSidebar";

interface NumericSummary {
  mean: number;
  median: number;
  min: number;
  max: number;
  [key: string]: number;
}

interface SelectionSummary {
  summary: Record<string, number> | NumericSummary;
  column_type: "numeric" | "categorical";
  total_selected: number;
}

interface UMAPProps {
  docState: DocState;
  setDocState: SetDocState;
  updateDocState: UpdateDocState;
  leftSidebarCollapsed: boolean;
}

export const UMAP: React.FC<UMAPProps> = ({
  docState,
  setDocState,
  leftSidebarCollapsed,
}) => {
  const containerRef = useRef<HTMLDivElement>(null);
  const svgRef = useRef<SVGSVGElement | null>(null);
  const umapData = docState.umap.data;
  const exportSVG = useExportSvg("umap", svgRef);
  const [loading, setLoading] = useState(false);
  const [dimensions, setDimensions] = useState({ width: 800, height: 600 });
  const [error, setError] = useState<string | null>(null);
  const [selectionSummary, setSelectionSummary] =
    useState<SelectionSummary | null>(null);
  const [selectedPoints, setSelectedPoints] = useState<string[]>([]);
  const [selectionActive, setSelectionActive] = useState(false);
  const [polygon, setPolygon] = useState<[number, number][]>([]);
  const loadingRef = useRef(false); // Prevent concurrent API calls

  useEffect(() => {
    exportSVG();
  }, [exportSVG]);

  // Use the reusable metadata colors hook
  const metadataColors = useMetadataColors({
    docId: docState.id,
    columnName: docState.umap.selectedMetadataColumn,
    enabled: docState.umap.colorBy === "metadata",
  });

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

  // Fetch UMAP data when parameters change (with debounce)
  useEffect(() => {
    if (!docState.id) return;

    const {
      n_neighbors,
      min_dist,
      minClusterSize,
      clusterEpsilon,
      params,
      data,
    } = docState.umap;

    // Check if parameters have changed
    const hasChanged =
      !params ||
      params.n_neighbors !== n_neighbors ||
      params.min_dist !== min_dist ||
      params.minClusterSize !== minClusterSize ||
      params.clusterEpsilon !== clusterEpsilon;

    // If data exists and parameters haven't changed, do nothing
    if (data && !hasChanged) {
      return;
    }

    // Add debounce to prevent multiple rapid API calls
    const timeoutId = setTimeout(async () => {
      // Skip if already loading
      if (loadingRef.current) {
        console.log("Skipping UMAP fetch - already loading");
        return;
      }

      loadingRef.current = true;
      setLoading(true);
      setError(null);

      try {
        const currentParams = {
          n_neighbors,
          min_dist,
          minClusterSize,
          clusterEpsilon,
        };

        const response = await window.pywebview.api.data.get_umap_data(
          docState.id,
          {
            n_neighbors,
            min_dist,
            threshold: minClusterSize,
            methods: [`hdbscan-${clusterEpsilon}`],
          },
        );

        setDocState((prev) => ({
          ...prev,
          umap: {
            ...prev.umap,
            data: response.data,
            params: currentParams,
          },
        }));
      } catch (err) {
        console.error("Error fetching UMAP data:", err);
        const errorMessage =
          err instanceof Error ? err.message : "Failed to load UMAP data";
        setError(errorMessage);
      } finally {
        setLoading(false);
        loadingRef.current = false;
      }
    }, 500); // 500ms debounce

    // Cleanup function to cancel pending requests
    return () => {
      clearTimeout(timeoutId);
    };
  }, [docState.id, docState.umap, setDocState]);

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

  // Main render effect - only renders the base visualization
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

        if (docState.umap.selectedMetadataColumn) {
          try {
            const summary =
              await window.pywebview.api.data.get_metadata_summary_for_selection(
                docState.id,
                selectedIds,
                docState.umap.selectedMetadataColumn,
              );
            setSelectionSummary(summary);
            console.log("Metadata summary for selection:", summary);
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
      .attr("r", 3) // Default size
      .attr("fill", "#1f77b4") // Default color
      .attr("opacity", 0.7) // Default opacity
      .attr("stroke", "white")
      .attr("stroke-width", 0.5)
      .attr("class", "umap-point") // Add class for easy selection
      .on("mouseover", (_event, d) => {
        tooltip.style("visibility", "visible");
        let content = `<strong>ID:</strong> ${d.id}`;
        if (d.cluster) {
          content += `<br><strong>Cluster:</strong> ${d.cluster}`;
        }
        if (docState.umap.colorBy === "metadata" && metadataColors.values) {
          const metadataValue = metadataColors.values[d.id];
          content += getMetadataTooltip(
            d.id,
            docState.umap.selectedMetadataColumn,
            metadataValue,
          );
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

    // Add brush
    const brush = d3
      .brush()
      .extent([
        [0, 0],
        [width, height],
      ])
      .filter((event) => event.ctrlKey && !event.button) // Only brush with Ctrl + Left-click
      .on("start brush", handleBrush)
      .on("end", handleBrushEnd);

    mainGroup
      .append("g")
      .attr("class", "brush")
      .attr("aria-label", "Selection brush area")
      .call(brush);

    // Prevent context menu on right-click to allow for other interactions
    svg.on("contextmenu", (event) => event.preventDefault());

    if (docState.umap.selectionMode === "polygon") {
      mainGroup.select(".brush").style("display", "none");
      svg.on("click", (event) => {
        const [mx, my] = d3.pointer(event);
        const newPolygon = [
          ...polygon,
          [mx - margin.left, my - margin.top],
        ] as [number, number][];
        setPolygon(newPolygon);
      });

      mainGroup
        .append("path")
        .datum(polygon)
        .attr("class", "selection-polygon")
        .attr("d", d3.line())
        .attr("stroke", "black")
        .attr("stroke-width", 2)
        .attr("fill", "rgba(0,0,0,0.1)");

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
          .on("click", async () => {
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

              // Get cluster summary
              const clusterCounts: Record<string, number> = {};
              for (const point of pointsInPolygon) {
                const cluster = point.cluster ?? 0;
                clusterCounts[cluster] = (clusterCounts[cluster] || 0) + 1;
              }

              // For metadata if available
              if (docState.umap.selectedMetadataColumn) {
                try {
                  const summary =
                    await window.pywebview.api.data.get_metadata_summary_for_selection(
                      docState.id,
                      selectedIds,
                      docState.umap.selectedMetadataColumn,
                    );
                  setSelectionSummary(summary);
                  console.log("Metadata summary for selection:", summary);
                } catch (err) {
                  console.error("Error fetching selection summary:", err);
                }
              }

              console.log(`Selected ${selectedIds.length} points`);
              console.log("Cluster distribution:", clusterCounts);
            }

            setPolygon([]);
          });
      }
    } else {
      mainGroup.select(".brush").style("display", "block");
      svg.on("click", null);
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
        xAxisGroup.call(xAxis.scale(currentXScale));
        yAxisGroup.call(yAxis.scale(currentYScale));

        // Update points positions
        svg
          .selectAll<SVGCircleElement, UMAPPoint>(".umap-point")
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

    // Cleanup on unmount
    return () => {
      tooltip.remove();
    };
  }, [
    umapData,
    dimensions,
    loading,
    error,
    docState.id,
    docState.umap.colorBy,
    docState.umap.selectedMetadataColumn,
    docState.umap.selectionMode,
    metadataColors.values,
    selectionActive,
    polygon,
  ]);

  // Update cluster stats panel separately
  useEffect(() => {
    if (!svgRef.current || !umapData) return;

    const svg = d3.select(svgRef.current);
    const mainGroup = svg.select("g");

    // Remove existing stats group
    mainGroup.select(".cluster-stats").remove();

    const clusterStats = umapData.clusterStats;
    const width = dimensions.width - 90; // margin.left + margin.right

    // Add cluster stats panel
    if (clusterStats && docState.umap.colorByCluster) {
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
  }, [umapData, dimensions, docState.umap.colorByCluster]);

  // Separate effect to update colors only
  useEffect(() => {
    if (!svgRef.current || !umapData) return;

    const svg = d3.select(svgRef.current);
    const points = svg.selectAll<SVGCircleElement, UMAPPoint>(".umap-point");

    if (points.empty()) return;

    // Get color function based on settings
    const getPointColor = (d: UMAPPoint): string => {
      if (docState.umap.colorBy === "metadata") {
        return metadataColors.getColor(d.id);
      }

      // Default cluster coloring
      if (!docState.umap.colorByCluster || !d.cluster) {
        return "#1f77b4"; // Default blue
      }
      const clusterId = d.cluster ?? 0;
      return clusterId === 0 ? "#cccccc" : distinctColor(clusterId); // Gray for noise points
    };

    // Update colors with a smooth transition
    points.transition().duration(300).attr("fill", getPointColor);

    // Update point size and opacity without transition
    points
      .attr("r", docState.umap.pointSize)
      .attr("opacity", docState.umap.opacity);

    // Restore selection state if there are selected points
    if (selectionActive && selectedPoints.length > 0) {
      points.classed("selected", (d) => selectedPoints.includes(d.id));
    }
  }, [
    docState.umap.colorBy,
    docState.umap.colorByCluster,
    docState.umap.pointSize,
    docState.umap.opacity,
    metadataColors.getColor,
    umapData,
    selectionActive,
    selectedPoints,
  ]);

  // Log selection info whenever it changes
  useEffect(() => {
    if (selectionSummary && selectedPoints.length > 0) {
      console.log(`Selected ${selectedPoints.length} points`);
      console.log("Metadata summary:", selectionSummary);
    }
  }, [selectionSummary, selectedPoints]);

  // Generate cluster distribution stats from selected points
  const getClusterDistribution = () => {
    if (!selectionActive || !selectedPoints.length || !umapData) return null;

    const clusterCounts: Record<string, number> = {};
    const filteredPoints = umapData.embedding.filter((point) =>
      selectedPoints.includes(point.id),
    );
    for (const point of filteredPoints) {
      const cluster = point.cluster ?? 0;
      clusterCounts[cluster] = (clusterCounts[cluster] || 0) + 1;
    }

    return clusterCounts;
  };

  const clusterDistribution = getClusterDistribution();

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
        {selectionActive && selectedPoints.length > 0 && (
          <div
            className="selection-info-panel"
            style={{
              position: "absolute",
              top: 10,
              right: 10,
              background: "rgba(255, 255, 255, 0.9)",
              padding: "10px",
              borderRadius: "4px",
              boxShadow:
                "0 1px 3px rgba(0,0,0,0.12), 0 1px 2px rgba(0,0,0,0.24)",
              zIndex: 100,
              maxWidth: "300px",
            }}
          >
            <h3 style={{ margin: "0 0 10px 0", fontSize: "16px" }}>
              Selection Summary
            </h3>
            <p style={{ margin: "5px 0" }}>
              <strong>Points selected:</strong> {selectedPoints.length}
            </p>

            {clusterDistribution && (
              <div>
                <h4 style={{ margin: "10px 0 5px 0", fontSize: "14px" }}>
                  Cluster Distribution:
                </h4>
                <ul style={{ margin: "0", paddingLeft: "20px" }}>
                  {Object.entries(clusterDistribution).map(
                    ([cluster, count]) => (
                      <li key={cluster} style={{ margin: "2px 0" }}>
                        <span
                          style={{
                            display: "inline-block",
                            width: "12px",
                            height: "12px",
                            backgroundColor:
                              cluster === "0"
                                ? "#cccccc"
                                : distinctColor(Number.parseInt(cluster)),
                            marginRight: "5px",
                            verticalAlign: "middle",
                          }}
                        />
                        Cluster {cluster === "0" ? "Noise" : cluster}: {count}{" "}
                        points
                      </li>
                    ),
                  )}
                </ul>
              </div>
            )}

            {selectionSummary && docState.umap.selectedMetadataColumn && (
              <div>
                <h4 style={{ margin: "10px 0 5px 0", fontSize: "14px" }}>
                  {docState.umap.selectedMetadataColumn} Summary:
                </h4>
                {selectionSummary.column_type === "numeric" ? (
                  <div>
                    <p style={{ margin: "2px 0" }}>
                      <strong>Mean:</strong>{" "}
                      {(
                        selectionSummary.summary as NumericSummary
                      ).mean?.toFixed(2)}
                    </p>
                    <p style={{ margin: "2px 0" }}>
                      <strong>Median:</strong>{" "}
                      {(
                        selectionSummary.summary as NumericSummary
                      ).median?.toFixed(2)}
                    </p>
                    <p style={{ margin: "2px 0" }}>
                      <strong>Min:</strong>{" "}
                      {(
                        selectionSummary.summary as NumericSummary
                      ).min?.toFixed(2)}
                    </p>
                    <p style={{ margin: "2px 0" }}>
                      <strong>Max:</strong>{" "}
                      {(
                        selectionSummary.summary as NumericSummary
                      ).max?.toFixed(2)}
                    </p>
                  </div>
                ) : (
                  <ul style={{ margin: "0", paddingLeft: "20px" }}>
                    {Object.entries(selectionSummary.summary).map(
                      ([value, count]) => (
                        <li key={value} style={{ margin: "2px 0" }}>
                          {value}: {count} points
                        </li>
                      ),
                    )}
                  </ul>
                )}
              </div>
            )}

            <button
              type="button"
              style={{
                marginTop: "10px",
                padding: "5px 10px",
                background: "#f44336",
                color: "white",
                border: "none",
                borderRadius: "4px",
                cursor: "pointer",
              }}
              onClick={() => {
                setSelectedPoints([]);
                setSelectionActive(false);
                setSelectionSummary(null);
                if (svgRef.current) {
                  d3.select(svgRef.current)
                    .selectAll(".umap-point")
                    .classed("selected", false);
                }
              }}
            >
              Clear Selection
            </button>
          </div>
        )}
        <svg
          ref={svgRef}
          id="umap-svg"
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
