import type React from "react";
import type { DocState } from "../../../appState";
import { distinctColor } from "../../../colors";

// Define NumericSummary interface locally as it's not exported from plotTypes
interface NumericSummary {
  mean: number;
  median: number;
  min: number;
  max: number;
  [key: string]: number;
}

interface SelectionSummaryProps {
  selectionActive: boolean;
  selectedPoints: string[];
  selectionSummary: {
    summary: Record<string, number> | NumericSummary;
    column_type: "numeric" | "categorical";
    total_selected: number;
  } | null;
  clusterDistribution: Record<string, number> | null;
  docState: DocState;
  updateSettings: (values: Partial<DocState["umap"]>) => void;
  onClearSelection: () => void;
  getMetadataColor: (id: string) => string;
  handleFileUpload: (
    event: React.ChangeEvent<HTMLInputElement>,
  ) => Promise<void>;
  uploading: boolean;
  uploadError: string | null;
  fileInputRef: React.RefObject<HTMLInputElement>;
}

export const SelectionSummaryPanel: React.FC<SelectionSummaryProps> = ({
  selectionActive,
  selectedPoints,
  selectionSummary,
  clusterDistribution,
  docState,
  updateSettings,
  onClearSelection,
  getMetadataColor,
  handleFileUpload,
  uploading,
  uploadError,
  fileInputRef,
}) => {
  // Function to handle color by selection change
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
    <div className="selection-info-panel">
      {selectionActive && selectedPoints.length > 0 ? (
        <>
          <h3>Selection Summary</h3>
          <p>
            <strong>Points selected:</strong> {selectedPoints.length}
          </p>

          <div style={{ marginBottom: "15px" }}>
            <select
              style={{
                width: "100%",
                padding: "5px",
                marginBottom: "10px",
                borderRadius: "4px",
                border: "1px solid #ccc",
              }}
              value={
                docState.umap.colorBy === "metadata"
                  ? docState.umap.selectedMetadataColumn || ""
                  : "cluster"
              }
              onChange={(e) => {
                const value = e.target.value;
                if (value === "cluster") {
                  updateSettings({
                    colorBy: "cluster",
                    colorByCluster: true,
                  });
                } else {
                  updateSettings({
                    colorBy: "metadata",
                    colorByCluster: false,
                    selectedMetadataColumn: value,
                  });
                }
              }}
            >
              <option value="cluster">Cluster</option>
              {docState.umap.uploadedMetadata?.columns.map((col) => (
                <option key={col} value={col}>
                  {col} ({docState.umap.uploadedMetadata?.columnTypes[col]})
                </option>
              ))}
            </select>
          </div>

          {/* Distribution visualization for current coloring method */}
          <div>
            <h4>Distribution:</h4>
            <div style={{ margin: "10px 0" }}>
              {docState.umap.colorBy === "cluster" &&
                clusterDistribution &&
                Object.entries(clusterDistribution).map(([cluster, count]) => {
                  const percentage = Math.round(
                    (count / selectedPoints.length) * 100,
                  );
                  const barWidth = `${percentage}%`;
                  return (
                    <div key={cluster} className="distribution-item">
                      <div className="color-indicator">
                        <span
                          className="color-dot"
                          style={{
                            backgroundColor:
                              cluster === "0"
                                ? "#cccccc"
                                : distinctColor(Number.parseInt(cluster)),
                          }}
                        />
                      </div>
                      <div className="item-label">
                        Cluster {cluster === "0" ? "Noise" : cluster}:
                      </div>
                      <div className="bar-container">
                        <div
                          className="bar"
                          style={{
                            width: barWidth,
                            backgroundColor:
                              cluster === "0"
                                ? "#cccccc"
                                : distinctColor(Number.parseInt(cluster)),
                          }}
                        />
                        <div className="bar-label">
                          {count} ({percentage}%)
                        </div>
                      </div>
                    </div>
                  );
                })}

              {docState.umap.colorBy === "metadata" &&
                selectionSummary &&
                (selectionSummary.column_type === "categorical" ? (
                  Object.entries(selectionSummary.summary).map(
                    ([value, count]) => {
                      const percentage = Math.round(
                        ((count as number) / selectedPoints.length) * 100,
                      );
                      const barWidth = `${percentage}%`;
                      return (
                        <div key={value} className="distribution-item">
                          <div className="item-name">{value}:</div>
                          <div className="bar-container">
                            <div
                              className="bar"
                              style={{
                                width: barWidth,
                                backgroundColor: getMetadataColor(value),
                              }}
                            />
                            <div className="bar-label">
                              {count as number} ({percentage}%)
                            </div>
                          </div>
                        </div>
                      );
                    },
                  )
                ) : (
                  <div>
                    <p>Numeric data summary:</p>
                    <ul className="summary-list">
                      <li>
                        Mean:{" "}
                        {(
                          selectionSummary.summary as NumericSummary
                        ).mean.toFixed(2)}
                      </li>
                      <li>
                        Median:{" "}
                        {(
                          selectionSummary.summary as NumericSummary
                        ).median.toFixed(2)}
                      </li>
                      <li>
                        Min:{" "}
                        {(
                          selectionSummary.summary as NumericSummary
                        ).min.toFixed(2)}
                      </li>
                      <li>
                        Max:{" "}
                        {(
                          selectionSummary.summary as NumericSummary
                        ).max.toFixed(2)}
                      </li>
                    </ul>
                  </div>
                ))}
            </div>
          </div>
          <button
            type="button"
            className="clear-selection-button"
            onClick={onClearSelection}
          >
            Clear Selection
          </button>
        </>
      ) : (
        <>
          {/* Section for metadata upload and selection controls */}
          <h3>UMAP Controls</h3>

          {/* Color Settings */}
          <div className="settings-section">
            <h4>Color Settings</h4>
            <label
              htmlFor="colorBySelect"
              style={{ display: "block", marginBottom: "5px" }}
            >
              Color Points By
            </label>
            <select
              id="colorBySelect"
              style={{
                width: "100%",
                padding: "5px",
                marginBottom: "10px",
                borderRadius: "4px",
                border: "1px solid #ccc",
              }}
              value={docState.umap.colorBy}
              onChange={(e) => {
                const value = e.target.value as "cluster" | "metadata";
                handleColorByChange(value);
              }}
            >
              <option value="cluster">Cluster</option>
              <option value="metadata">Metadata</option>
            </select>

            {docState.umap.colorBy === "metadata" && (
              <>
                <div style={{ marginTop: "10px" }}>
                  <label
                    htmlFor="metadata-upload"
                    style={{ display: "block", marginBottom: "5px" }}
                  >
                    Upload Metadata CSV
                  </label>
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

                {docState.umap.uploadedMetadata && (
                  <>
                    <div style={{ marginTop: "10px" }}>
                      <label
                        htmlFor="metadataColumnSelect"
                        style={{ display: "block", marginBottom: "5px" }}
                      >
                        Metadata Column
                      </label>
                      <select
                        id="metadataColumnSelect"
                        style={{
                          width: "100%",
                          padding: "5px",
                          borderRadius: "4px",
                          border: "1px solid #ccc",
                        }}
                        value={docState.umap.selectedMetadataColumn || ""}
                        onChange={(e) =>
                          handleMetadataColumnChange(e.target.value)
                        }
                      >
                        {docState.umap.uploadedMetadata.columns.map((col) => (
                          <option key={col} value={col}>
                            {col} (
                            {docState.umap.uploadedMetadata?.columnTypes[col]})
                          </option>
                        ))}
                      </select>
                    </div>

                    <div className="metadata-stats">
                      <h5>Match Statistics</h5>
                      <p>
                        Matched:{" "}
                        {docState.umap.uploadedMetadata.matchStats
                          .exactMatches +
                          docState.umap.uploadedMetadata.matchStats
                            .versionMatches}{" "}
                        /{" "}
                        {
                          docState.umap.uploadedMetadata.matchStats
                            .totalSequenceIds
                        }{" "}
                        (
                        {
                          docState.umap.uploadedMetadata.matchStats
                            .matchPercentage
                        }
                        %)
                      </p>
                      <p className="match-details">
                        Exact:{" "}
                        {docState.umap.uploadedMetadata.matchStats.exactMatches}
                      </p>
                      <p className="match-details">
                        Version:{" "}
                        {
                          docState.umap.uploadedMetadata.matchStats
                            .versionMatches
                        }
                      </p>
                      <p className="match-details">
                        Unmatched:{" "}
                        {docState.umap.uploadedMetadata.matchStats.unmatched}
                      </p>
                    </div>
                  </>
                )}
              </>
            )}
          </div>

          {/* Visual Settings */}
          <div className="settings-section">
            <h4>Visual Settings</h4>

            <div style={{ marginBottom: "10px" }}>
              <label
                htmlFor="pointSize"
                style={{ display: "block", marginBottom: "5px" }}
              >
                Point Size: {docState.umap.pointSize}
              </label>
              <input
                type="range"
                id="pointSize"
                min="1"
                max="20"
                step="0.5"
                value={docState.umap.pointSize}
                onChange={(e) =>
                  updateSettings({
                    pointSize: Number.parseFloat(e.target.value),
                  })
                }
                style={{ width: "100%" }}
              />
            </div>

            <div>
              <label
                htmlFor="opacity"
                style={{ display: "block", marginBottom: "5px" }}
              >
                Opacity: {docState.umap.opacity.toFixed(2)}
              </label>
              <input
                type="range"
                id="opacity"
                min="0.1"
                max="1"
                step="0.05"
                value={docState.umap.opacity}
                onChange={(e) =>
                  updateSettings({
                    opacity: Number.parseFloat(e.target.value),
                  })
                }
                style={{ width: "100%" }}
              />
            </div>
          </div>

          {/* Help Button */}
          <button type="button" className="tooltip-container help-button">
            <span className="help-icon">ⓘ</span>
            <span>Usage Instructions</span>
            <div className="tooltip-content">
              <p>
                <strong>Navigation:</strong> Left-click drag to pan, scroll to
                zoom
              </p>
              <p>
                <strong>Brush Selection:</strong> Right-click drag to select
                points with a rectangular brush
              </p>
              <p>
                <strong>Polygon Selection:</strong> Hold Ctrl + right-click to
                place polygon vertices, then Ctrl + right-click on the red
                starting point to complete selection
              </p>
              <p>
                <strong>Note:</strong> Both selection methods are always
                available
              </p>
            </div>
          </button>
        </>
      )}
    </div>
  );
};
