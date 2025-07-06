import React from "react";
import type { MetaData } from "../plotTypes";

export const AlignmentStats = ({
  metaData,
  dataLength,
  activeDataSet,
}: {
  metaData?: MetaData;
  dataLength?: number;
  activeDataSet: string;
}) => {
  const [isExpanded, setIsExpanded] = React.useState(false);

  const isLzani = metaData?.run?.analysis_method === "lzani";
  const unalignedCount = metaData?.unaligned_count || 0;
  const alignedCount = dataLength || 0;
  const totalCount = alignedCount + unalignedCount;

  const alignedPercent = totalCount > 0 ? (alignedCount / totalCount) * 100 : 0;
  const unalignedPercent =
    totalCount > 0 ? (unalignedCount / totalCount) * 100 : 0;

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
  if (!stats) {
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

          {/* Mini bar chart - only for LZANI scores */}
          {activeDataSet === "scores" && isLzani && unalignedCount > 0 && (
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

          {/* Distribution Statistics */}
          {stats && (
            <div
              style={{
                marginTop:
                  activeDataSet === "scores" && isLzani && unalignedCount > 0
                    ? "15px"
                    : "0",
                borderTop:
                  activeDataSet === "scores" && isLzani && unalignedCount > 0
                    ? "1px solid #ddd"
                    : "none",
                paddingTop:
                  activeDataSet === "scores" && isLzani && unalignedCount > 0
                    ? "10px"
                    : "0",
              }}
            >
              <div
                style={{
                  fontSize: "13px",
                  fontWeight: "bold",
                  marginBottom: "8px",
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
        </div>
      )}
    </>
  );
};
