import React from "react";
import type { ClusterStats, MetaData } from "../../plotTypes";

type ClusterStatsProps = {
  metaData?: MetaData;
  clusterStats?: ClusterStats;
};

export const ClusterStatsDisplay = ({
  metaData,
  clusterStats,
}: ClusterStatsProps) => {
  const [isExpanded, setIsExpanded] = React.useState(false);

  const isLzani = metaData?.run?.analysis_method === "lzani";
  const unalignedCount = metaData?.unaligned_count || 0;
  const alignedCount = metaData?.distribution_stats?.count || 0;
  const totalCount = alignedCount + unalignedCount;

  const alignedPercent = totalCount > 0 ? (alignedCount / totalCount) * 100 : 0;
  const unalignedPercent =
    totalCount > 0 ? (unalignedCount / totalCount) * 100 : 0;

  if (!clusterStats) {
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
          title="Show cluster statistics"
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
              Cluster Summary
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

          {/* Mini bar chart - only for LZANI scores */}
          {isLzani && unalignedCount > 0 && (
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
                      Unaligned: {unalignedCount.toLocaleString()} (
                      {unalignedPercent.toFixed(1)}%)
                    </span>
                  </div>
                </div>
              </div>

              <div
                style={{ fontSize: "11px", color: "#999", marginTop: "10px" }}
              >
                Pairs below 70% threshold are excluded
              </div>
            </>
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
