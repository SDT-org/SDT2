import React from "react";
import type {
  ClusterStats,
  HeatmapData,
  HeatmapSettings,
  MetaData,
} from "../plotTypes";

type AlignmentStatsProps = {
  metaData?: MetaData;
  dataLength?: number;
  activeDataSet: string;
  clusterStats?: ClusterStats;
  heatmapSettings?: HeatmapSettings;
  rawData?: HeatmapData;
};

export const AlignmentStats = ({
  metaData,
  dataLength,
  activeDataSet,
  clusterStats,
  heatmapSettings,
  rawData,
}: AlignmentStatsProps) => {
  const [isExpanded, setIsExpanded] = React.useState(false);

  const isLzani = metaData?.run?.analysis_method === "lzani";
  const isMaskingEnabled = heatmapSettings?.hideValuesBelowEnabled || false;
  const maskThreshold = heatmapSettings?.hideValuesBelow || 0;

  // Calculate masked values count if masking is enabled
  let maskedCount = 0;
  if (isMaskingEnabled && rawData && maskThreshold > 0) {
    for (let i = 0; i < rawData.length; i++) {
      for (let j = 0; j < i; j++) {
        const value = rawData[i]?.[j];
        if (
          value !== null &&
          value !== undefined &&
          value !== 0 &&
          value < maskThreshold
        ) {
          maskedCount++;
        }
      }
    }
  }

  // For LZANI, unaligned count comes from backend
  // For both LZANI and Parasail, add masked values to unaligned count
  const baseUnalignedCount = metaData?.unaligned_count || 0;
  const totalUnalignedCount = baseUnalignedCount + maskedCount;

  // Aligned count should exclude masked values
  const alignedCount = (dataLength || 0) - maskedCount;
  const totalCount = alignedCount + totalUnalignedCount;

  const alignedPercent = totalCount > 0 ? (alignedCount / totalCount) * 100 : 0;
  const unalignedPercent =
    totalCount > 0 ? (totalUnalignedCount / totalCount) * 100 : 0;

  // Show stats if:
  // 1. LZANI with unaligned values
  // 2. Any method with masking enabled
  const showAlignmentStats =
    (isLzani && baseUnalignedCount > 0) ||
    (isMaskingEnabled && maskedCount > 0);

  // Get the appropriate stats based on active dataset
  const stats =
    activeDataSet === "scores"
      ? metaData?.distribution_stats
      : activeDataSet === "gc"
        ? metaData?.gc_stats
        : activeDataSet === "length"
          ? metaData?.length_stats
          : null;

  // Don't show if no stats are available
  if (!stats && !clusterStats) {
    return null;
  }

  return (
    <>
      {/* Collapsed button */}
      {!isExpanded && (
        <button
          type="button"
          onClick={() => setIsExpanded(true)}
          style={{
            position: "absolute",
            bottom: "20px",
            right: "20px",
            width: "40px",
            height: "40px",
            backgroundColor: "#f8f8f8",
            border: "1px solid #ddd",
            borderRadius: "4px",
            cursor: "pointer",
            display: "flex",
            alignItems: "center",
            justifyContent: "center",
            fontSize: "18px",
            color: "#666",
            boxShadow: "0 2px 4px rgba(0,0,0,0.1)",
            transition: "all 0.2s ease",
            zIndex: 1000,
          }}
          onMouseEnter={(e) => {
            e.currentTarget.style.backgroundColor = "#e8e8e8";
          }}
          onMouseLeave={(e) => {
            e.currentTarget.style.backgroundColor = "#f8f8f8";
          }}
          title="Show alignment statistics"
        >
          📊
        </button>
      )}

      {/* Expanded panel */}
      {isExpanded && (
        <div
          style={{
            position: "absolute",
            bottom: "20px",
            right: "20px",
            width: "250px",
            backgroundColor: "#f8f8f8",
            border: "1px solid #ddd",
            borderRadius: "4px",
            padding: "15px",
            boxShadow: "0 2px 8px rgba(0,0,0,0.15)",
            zIndex: 1000,
          }}
        >
          <div
            style={{
              display: "flex",
              justifyContent: "space-between",
              alignItems: "center",
              marginBottom: "10px",
            }}
          >
            <div style={{ fontSize: "14px", fontWeight: "bold", margin: 0 }}>
              {activeDataSet === "scores"
                ? "Alignment Summary"
                : activeDataSet === "gc"
                  ? "GC Content Summary"
                  : "Sequence Length Summary"}
            </div>
            <button
              type="button"
              onClick={() => setIsExpanded(false)}
              style={{
                background: "none",
                border: "none",
                cursor: "pointer",
                fontSize: "16px",
                color: "#666",
                padding: "0 4px",
              }}
              title="Close"
            >
              ✕
            </button>
          </div>

          {/* Mini bar chart - show for LZANI or when masking is enabled */}
          {activeDataSet === "scores" && showAlignmentStats && (
            <>
              <div style={{ marginBottom: "10px" }}>
                <div
                  style={{
                    display: "flex",
                    height: "30px",
                    marginBottom: "5px",
                  }}
                >
                  <div
                    style={{
                      width: `${alignedPercent}%`,
                      backgroundColor: "#4CAF50",
                      display: "flex",
                      alignItems: "center",
                      justifyContent: "center",
                      color: "white",
                      fontSize: "12px",
                      minWidth: alignedPercent > 10 ? "auto" : "0",
                    }}
                  >
                    {alignedPercent > 10 && `${alignedPercent.toFixed(1)}%`}
                  </div>
                  <div
                    style={{
                      width: `${unalignedPercent}%`,
                      backgroundColor: "#f44336",
                      display: "flex",
                      alignItems: "center",
                      justifyContent: "center",
                      color: "white",
                      fontSize: "12px",
                      minWidth: unalignedPercent > 10 ? "auto" : "0",
                    }}
                  >
                    {unalignedPercent > 10 && `${unalignedPercent.toFixed(1)}%`}
                  </div>
                </div>

                {/* Legend */}
                <div style={{ fontSize: "12px", color: "#666" }}>
                  <div
                    style={{
                      display: "flex",
                      alignItems: "center",
                      marginBottom: "3px",
                    }}
                  >
                    <div
                      style={{
                        width: "12px",
                        height: "12px",
                        backgroundColor: "#4CAF50",
                        marginRight: "5px",
                      }}
                    />
                    <span>
                      Aligned: {alignedCount.toLocaleString()} (
                      {alignedPercent.toFixed(1)}%)
                    </span>
                  </div>
                  <div style={{ display: "flex", alignItems: "center" }}>
                    <div
                      style={{
                        width: "12px",
                        height: "12px",
                        backgroundColor: "#f44336",
                        marginRight: "5px",
                      }}
                    />
                    <span>
                      Unaligned: {totalUnalignedCount.toLocaleString()} (
                      {unalignedPercent.toFixed(1)}%)
                    </span>
                  </div>
                </div>
              </div>

              <div
                style={{ fontSize: "11px", color: "#999", marginTop: "10px" }}
              >
                {isLzani && baseUnalignedCount > 0 && maskedCount > 0
                  ? `Pairs below 70% (LZANI) and ${maskThreshold}% (masked) are excluded`
                  : isLzani && baseUnalignedCount > 0
                    ? "Pairs below 70% threshold are excluded"
                    : `Pairs below ${maskThreshold}% threshold are excluded`}
              </div>
            </>
          )}

          {/* Distribution Statistics */}
          {stats && (
            <div
              style={{
                marginTop:
                  activeDataSet === "scores" && showAlignmentStats
                    ? "15px"
                    : "0",
                borderTop:
                  activeDataSet === "scores" && showAlignmentStats
                    ? "1px solid #ddd"
                    : "none",
                paddingTop:
                  activeDataSet === "scores" && showAlignmentStats
                    ? "10px"
                    : "0",
              }}
            >
              <div
                style={{
                  fontSize: "13px",
                  fontWeight: "bold",
                  marginBottom: "8px",
                  color: activeDataSet === "scores" ? "#d48b91" : "#333",
                }}
              >
                Distribution Statistics
              </div>
              <div
                style={{ fontSize: "11px", color: "#666", lineHeight: "1.4" }}
              >
                {(() => {
                  const suffix = activeDataSet === "length" ? " nt" : "%";
                  const decimals = activeDataSet === "length" ? 0 : 2;

                  return (
                    <>
                      <div>
                        Mean: {stats.mean.toFixed(decimals)}
                        {suffix}
                      </div>
                      <div>
                        Median: {stats.median.toFixed(decimals)}
                        {suffix}
                      </div>
                      <div>
                        Std Dev: {stats.std.toFixed(decimals)}
                        {activeDataSet === "length" ? " nt" : ""}
                      </div>
                      <div>
                        Min: {stats.min.toFixed(decimals)}
                        {suffix}
                      </div>
                      <div>
                        Max: {stats.max.toFixed(decimals)}
                        {suffix}
                      </div>
                      <div>
                        Q1: {stats.q1.toFixed(decimals)}
                        {suffix}
                      </div>
                      <div>
                        Q3: {stats.q3.toFixed(decimals)}
                        {suffix}
                      </div>
                    </>
                  );
                })()}
              </div>
            </div>
          )}

          {/* Cluster Statistics */}
          {clusterStats && (
            <div
              style={{
                marginTop: "15px",
                borderTop: "1px solid #ddd",
                paddingTop: "10px",
              }}
            >
              <div
                style={{
                  fontSize: "13px",
                  fontWeight: "bold",
                  marginBottom: "8px",
                  color: "#333",
                }}
              >
                Cluster Statistics
              </div>
              <div
                style={{ fontSize: "11px", color: "#666", lineHeight: "1.4" }}
              >
                <div>Total Clusters: {clusterStats.total_clusters}</div>
                <div>Largest Cluster: {clusterStats.largest_cluster}</div>
                <div>Singleton Clusters: {clusterStats.singleton_clusters}</div>
              </div>
            </div>
          )}
        </div>
      )}
    </>
  );
};
