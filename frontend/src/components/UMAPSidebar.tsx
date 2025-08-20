import type React from "react";
import { useRef, useState } from "react";
import { useAppState } from "../appState";
import { useDocState } from "../hooks/useDocState";
import type { UMAPSettings } from "../plotTypes";
import { Select, SelectItem } from "./primitives/Select";
import { Slider } from "./primitives/Slider";

interface UMAPSidebarProps {
  settings: UMAPSettings;
  updateSettings: (values: Partial<UMAPSettings>) => void;
  leftSidebarCollapsed: boolean;
}

export const UMAPSidebar: React.FC<UMAPSidebarProps> = ({
  settings,
  updateSettings,
}) => {
  const { appState, setAppState } = useAppState();
  const { docState } = useDocState(
    appState.activeDocumentId,
    appState,
    setAppState,
  );
  const [uploading, setUploading] = useState(false);
  const [uploadError, setUploadError] = useState<string | null>(null);
  const fileInputRef = useRef<HTMLInputElement>(null);

  const handleParameterChange = (param: string, value: number | boolean) => {
    updateSettings({ [param]: value });
  };

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

  // Handle color by selection change
  const handleColorByChange = (value: string) => {
    updateSettings({
      colorBy: value as "cluster" | "metadata",
      colorByCluster: value === "cluster", // Keep backward compatibility
    });
  };

  // Handle metadata column selection
  const handleMetadataColumnChange = (value: string) => {
    updateSettings({ selectedMetadataColumn: value });
  };

  return (
    <div className="app-sidebar app-sidebar-right umap-sidebar">
      <div className="app-sidebar-body">
        <div className="sidebar-section">
          <h3>UMAP Parameters</h3>
          <p className="sidebar-description">
            UMAP (Uniform Manifold Approximation and Projection) reduces
            high-dimensional distance data to 2D for visualization.
          </p>

          <div className="sidebar-item">
            <label htmlFor="n_neighbors">n_neighbors</label>
            <p className="param-description">
              Controls how UMAP balances local vs global structure. Lower values
              focus on local structure.
            </p>
            <Slider
              id="n_neighbors"
              value={settings.n_neighbors}
              onChangeEnd={(value: number) =>
                handleParameterChange("n_neighbors", value)
              }
              minValue={2}
              maxValue={200}
              step={1}
            />
            <span className="value-display">{settings.n_neighbors}</span>
          </div>

          <div className="sidebar-item">
            <label htmlFor="min_dist">min_dist</label>
            <p className="param-description">
              Minimum distance between points. Lower values create tighter
              clusters.
            </p>
            <Slider
              id="min_dist"
              value={settings.min_dist}
              onChangeEnd={(value: number) =>
                handleParameterChange("min_dist", value)
              }
              minValue={0}
              maxValue={1}
              step={0.01}
            />
            <span className="value-display">
              {settings.min_dist.toFixed(2)}
            </span>
          </div>
        </div>

        <div className="sidebar-section">
          <h3>HDBSCAN Clustering</h3>
          <p className="sidebar-description">
            HDBSCAN finds clusters of varying densities in the data. Points not
            belonging to any cluster are marked as noise.
          </p>

          <div className="sidebar-item">
            <label htmlFor="minClusterSize">Min Cluster Size</label>
            <p className="param-description">
              Minimum number of points required to form a cluster. Smaller
              values find more clusters.
            </p>
            <Slider
              id="minClusterSize"
              value={settings.minClusterSize}
              onChangeEnd={(value: number) =>
                handleParameterChange("minClusterSize", value)
              }
              minValue={2}
              maxValue={50}
              step={1}
            />
            <span className="value-display">{settings.minClusterSize}</span>
          </div>

          <div className="sidebar-item">
            <label htmlFor="clusterEpsilon">
              Cluster Selection Epsilon (%)
            </label>
            <p className="param-description">
              Minimum similarity percentage for clustering. Points must be at
              least this similar to belong to the same cluster. 0 uses the full
              hierarchy.
            </p>
            <Slider
              id="clusterEpsilon"
              value={settings.clusterEpsilon}
              onChangeEnd={(value: number) =>
                handleParameterChange("clusterEpsilon", value)
              }
              minValue={0}
              maxValue={100}
              step={0.5}
            />
            <span className="value-display">
              {settings.clusterEpsilon.toFixed(1)}
            </span>
          </div>

          <div className="sidebar-item">
            <Select
              id="colorBy"
              label="Color Points By"
              selectedKey={settings.colorBy}
              onSelectionChange={(key) => handleColorByChange(key as string)}
            >
              <SelectItem id="cluster">Cluster</SelectItem>
              <SelectItem id="metadata">Metadata</SelectItem>
            </Select>
          </div>

          {settings.colorBy === "metadata" && (
            <>
              <div className="sidebar-item">
                <label htmlFor="metadata-upload">Upload Metadata CSV</label>
                <input
                  id="metadata-upload"
                  ref={fileInputRef}
                  type="file"
                  accept=".csv"
                  onChange={handleFileUpload}
                  disabled={uploading}
                  style={{ width: "100%" }}
                />
                {uploading && <p className="upload-status">Uploading...</p>}
                {uploadError && <p className="upload-error">{uploadError}</p>}
              </div>

              {settings.uploadedMetadata && (
                <>
                  <div className="sidebar-item">
                    <Select
                      id="metadataColumn"
                      label="Metadata Column"
                      selectedKey={
                        settings.selectedMetadataColumn ||
                        settings.uploadedMetadata.columns[0] ||
                        null
                      }
                      onSelectionChange={(key) =>
                        key && handleMetadataColumnChange(key as string)
                      }
                    >
                      {settings.uploadedMetadata.columns.map((col) => (
                        <SelectItem key={col} id={col}>
                          {col} ({settings.uploadedMetadata?.columnTypes[col]})
                        </SelectItem>
                      ))}
                    </Select>
                  </div>

                  <div className="metadata-stats">
                    <h4>Match Statistics</h4>
                    <p>
                      Matched:{" "}
                      {settings.uploadedMetadata.matchStats.exactMatches +
                        settings.uploadedMetadata.matchStats
                          .versionMatches}{" "}
                      / {settings.uploadedMetadata.matchStats.totalSequenceIds}{" "}
                      ( {settings.uploadedMetadata.matchStats.matchPercentage}%)
                    </p>
                    <p className="match-details">
                      Exact: {settings.uploadedMetadata.matchStats.exactMatches}
                    </p>
                    <p className="match-details">
                      Version:{" "}
                      {settings.uploadedMetadata.matchStats.versionMatches}
                    </p>
                    <p className="match-details">
                      Unmatched:{" "}
                      {settings.uploadedMetadata.matchStats.unmatched}
                    </p>
                  </div>
                </>
              )}
            </>
          )}
        </div>

        <div className="sidebar-section">
          <h3>Visual Settings</h3>

          <div className="sidebar-item">
            <Select
              id="selectionMode"
              label="Selection Mode"
              selectedKey={settings.selectionMode}
              onSelectionChange={(key) =>
                updateSettings({ selectionMode: key as "brush" | "polygon" })
              }
            >
              <SelectItem id="brush">Brush</SelectItem>
              <SelectItem id="polygon">Polygon</SelectItem>
            </Select>
          </div>

          <div className="sidebar-item">
            <label htmlFor="pointSize">Point Size</label>
            <Slider
              id="pointSize"
              value={settings.pointSize}
              onChangeEnd={(value: number) =>
                handleParameterChange("pointSize", value)
              }
              minValue={1}
              maxValue={20}
              step={0.5}
            />
            <span className="value-display">{settings.pointSize}</span>
          </div>

          <div className="sidebar-item">
            <label htmlFor="opacity">Opacity</label>
            <Slider
              id="opacity"
              value={settings.opacity}
              onChangeEnd={(value: number) =>
                handleParameterChange("opacity", value)
              }
              minValue={0.1}
              maxValue={1}
              step={0.05}
            />
            <span className="value-display">{settings.opacity.toFixed(2)}</span>
          </div>
        </div>
      </div>
    </div>
  );
};
