import * as d3 from "d3";
import React, { useEffect, useRef, useState } from "react";
import type { DocState, SetDocState, UpdateDocState } from "../../../appState";
import { useExportSvg } from "../../../hooks/useExportSvg";
import { useMetadataColors } from "../../../hooks/useMetadataColors";
import { D3SvgUMAP } from "./D3SvgUMAP";
import { SelectionSummaryPanel } from "./SelectionSummaryPanel";
import { UMAPSidebar } from "./UMAPSidebar";
import "./umap.css";

interface SelectionSummary {
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
  const fileInputRef = useRef<HTMLInputElement>(null);
  const umapData = docState.umap.data;
  const exportSVG = useExportSvg("umap", svgRef);
  const [loading, setLoading] = useState(false);
  const [uploading, setUploading] = useState(false);
  const [uploadError, setUploadError] = useState<string | null>(null);
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

  // Handle file upload for metadata
  const handleFileUpload = async (
    event: React.ChangeEvent<HTMLInputElement>,
  ) => {
    const file = event.target.files?.[0];
    if (!file) return;

    setUploading(true);
    setUploadError(null);

    try {
      const content = await file.text();
      const response = await window.pywebview.api.data.upload_metadata(
        docState.id,
        content,
      );

      updateSettings({
        uploadedMetadata: {
          columns: response.columns,
          columnTypes: response.column_types,
          matchStats: {
            totalMetadataIds: response.match_stats.total_metadata_ids,
            totalSequenceIds: response.match_stats.total_sequence_ids,
            exactMatches: response.match_stats.exact_matches,
            versionMatches: response.match_stats.version_matches,
            unmatched: response.match_stats.unmatched,
            matchPercentage: response.match_stats.match_percentage,
          },
        },
        colorBy: "metadata",
        ...(response.columns.length > 0 && {
          selectedMetadataColumn: response.columns[0],
        }),
      });
    } catch (error) {
      setUploadError(error instanceof Error ? error.message : "Upload failed");
    } finally {
      setUploading(false);
      if (fileInputRef.current) {
        fileInputRef.current.value = "";
      }
    }
  };

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

  // Update selection summary when metadata column changes
  useEffect(() => {
    // Only fetch new summary if we have selected points and a metadata column is selected
    if (
      selectedPoints.length > 0 &&
      docState.umap.selectedMetadataColumn &&
      selectionActive
    ) {
      const fetchSummary = async () => {
        try {
          // Ensure the column name is not undefined before passing to API
          const columnName = docState.umap.selectedMetadataColumn;
          if (!columnName) return;

          const summary =
            await window.pywebview.api.data.get_metadata_summary_for_selection(
              docState.id,
              selectedPoints,
              columnName,
            );
          setSelectionSummary(summary);
          console.log("Updated metadata summary for selection:", summary);
        } catch (err) {
          console.error("Error fetching updated selection summary:", err);
        }
      };

      fetchSummary();
    }
  }, [
    docState.id,
    docState.umap.selectedMetadataColumn,
    selectedPoints,
    selectionActive,
  ]);

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

  // Handle clearing selection
  const handleClearSelection = () => {
    setSelectionActive(false);
    setSelectionSummary(null);
    setSelectedPoints([]);
    const svg = svgRef.current;
    if (svg) {
      const d3svg = d3.select(svg);
      d3svg.selectAll(".umap-point").classed("selected", false);
    }
  };

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
        {error && <div className="error-message">{error}</div>}

        {/* Selection Summary Panel Component */}
        <SelectionSummaryPanel
          selectionActive={selectionActive}
          selectedPoints={selectedPoints}
          selectionSummary={selectionSummary}
          clusterDistribution={clusterDistribution}
          docState={docState}
          updateSettings={updateSettings}
          onClearSelection={handleClearSelection}
          getMetadataColor={metadataColors.getColor}
          handleFileUpload={handleFileUpload}
          uploading={uploading}
          uploadError={uploadError}
          fileInputRef={fileInputRef}
        />

        {/* D3 UMAP Visualization Component */}
        {umapData && (
          <D3SvgUMAP
            umapData={umapData}
            dimensions={dimensions}
            pointSize={docState.umap.pointSize}
            opacity={docState.umap.opacity}
            colorBy={docState.umap.colorBy}
            colorByCluster={docState.umap.colorByCluster}
            selectedMetadataColumn={docState.umap.selectedMetadataColumn || ""}
            metadataValues={metadataColors.values || {}}
            selectedPoints={selectedPoints}
            selectionActive={selectionActive}
            polygon={polygon}
            svgRef={svgRef}
            getMetadataColor={metadataColors.getColor}
            setPolygon={setPolygon}
            setSelectedPoints={setSelectedPoints}
            setSelectionActive={setSelectionActive}
            setSelectionSummary={setSelectionSummary}
            docId={docState.id}
            error={error}
            loading={loading}
          />
        )}
      </div>

      {/* UMAP Sidebar Component */}
      <UMAPSidebar
        settings={docState.umap}
        updateSettings={updateSettings}
        leftSidebarCollapsed={leftSidebarCollapsed}
      />
    </>
  );
};
