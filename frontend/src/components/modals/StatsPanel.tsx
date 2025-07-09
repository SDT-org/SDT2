import React from "react";
import type { ClusterDataItem, ClusterStats, MetaData } from "../../plotTypes";

type StatsPanelProps = {
  metaData?: MetaData;
  dataLength?: number;
  activeDataSet?: string;
  clusterStats?: ClusterStats;
  clusterData?: ClusterDataItem[] | undefined;
  title?: string;
  panelType: "alignment" | "cluster" | "distribution";
};

export const StatsPanel = ({
  metaData,
  dataLength,
  activeDataSet = "scores",
  clusterStats,
  clusterData,
  title,
  panelType,
}: StatsPanelProps) => {
  const [isExpanded, setIsExpanded] = React.useState(false);
  const [isBarCollapsed, setIsBarCollapsed] = React.useState(false);

  const isLzani = metaData?.run?.analysis_method === "lzani";
  const unalignedCount = metaData?.unaligned_count || 0;
  const alignedCount = dataLength || metaData?.distribution_stats?.count || 0;
  const totalCount = alignedCount + unalignedCount;

  const alignedPercent = totalCount > 0 ? (alignedCount / totalCount) * 100 : 0;
  const unalignedPercent =
    totalCount > 0 ? (unalignedCount / totalCount) * 100 : 0;

  // Get panel type to match active dataset
  const getStats = () => {
    if (panelType === "distribution" || panelType === "alignment") {
      return activeDataSet === "scores"
        ? metaData?.distribution_stats
        : activeDataSet === "gc"
          ? metaData?.gc_stats
          : activeDataSet === "length"
            ? metaData?.length_stats
            : null;
    }
    return null;
  };

  const stats = getStats();

  //only lzani has unaligned fraction
  const shouldShowPanel =
    stats || clusterStats || (isLzani && unalignedCount > 0);
  if (!shouldShowPanel) {
    return null;
  }

  const getPanelTitle = () => {
    if (title) return title;

    if (panelType === "cluster") {
      return "Cluster Summary";
    }
    if (panelType === "distribution" || panelType === "alignment") {
      return activeDataSet === "scores"
        ? "Alignment Summary"
        : activeDataSet === "gc"
          ? "GC Content Summary"
          : "Sequence Length Summary";
    }

    return "Statistics";
  };

  const renderAlignmentBar = () => {
    if (!isLzani || unalignedCount === 0 || activeDataSet !== "scores") {
      return null;
    }

    return (
      <div className="stats-panel-alignment">
        <div className="stats-panel-alignment-header">
          <span>Alignment Status</span>
          <button
            type="button"
            onClick={() => setIsBarCollapsed(!isBarCollapsed)}
            title={isBarCollapsed ? "Expand" : "Collapse"}
          >
            {isBarCollapsed ? "▶" : "▼"}
          </button>
        </div>

        {!isBarCollapsed && (
          <>
            <div className="stats-panel-alignment-bar">
              <div
                className="stats-panel-alignment-bar-aligned"
                style={{
                  width: `${alignedPercent}%`,
                  minWidth: alignedPercent > 10 ? "auto" : "0",
                }}
              >
                {alignedPercent > 10 && `${alignedPercent.toFixed(1)}%`}
              </div>
              <div
                className="stats-panel-alignment-bar-unaligned"
                style={{
                  width: `${unalignedPercent}%`,
                  minWidth: unalignedPercent > 10 ? "auto" : "0",
                }}
              >
                {unalignedPercent > 10 && `${unalignedPercent.toFixed(1)}%`}
              </div>
            </div>

            {/* Legend */}
            <div className="stats-panel-alignment-legend">
              <div className="stats-panel-alignment-legend-item">
                <div className="stats-panel-alignment-legend-item-color aligned" />
                <span>
                  Aligned: {alignedCount.toLocaleString()} (
                  {alignedPercent.toFixed(1)}%)
                </span>
              </div>
              <div className="stats-panel-alignment-legend-item">
                <div className="stats-panel-alignment-legend-item-color unaligned" />
                <span>
                  Unaligned: {unalignedCount.toLocaleString()} (
                  {unalignedPercent.toFixed(1)}%)
                </span>
              </div>
            </div>
          </>
        )}
      </div>
    );
  };

  const renderDistributionStats = () => {
    if (!stats) return null;

    const suffix = activeDataSet === "length" ? " nt" : "%";
    const decimals = activeDataSet === "length" ? 0 : 2;

    return (
      <div className="stats-panel-section">
        <div className="stats-panel-section-title">Distribution</div>
        <div className="stats-panel-section-content">
          <div className="stats-panel-stat-row">
            <span>Mean:</span>
            <span>
              {stats.mean.toFixed(decimals)}
              {suffix}
            </span>
          </div>
          <div className="stats-panel-stat-row">
            <span>Median:</span>
            <span>
              {stats.median.toFixed(decimals)}
              {suffix}
            </span>
          </div>
          <div className="stats-panel-stat-row">
            <span>Std Dev:</span>
            <span>
              {stats.std.toFixed(decimals)}
              {activeDataSet === "length" ? " nt" : ""}
            </span>
          </div>
          <div className="stats-panel-stat-row">
            <span>Min:</span>
            <span>
              {stats.min.toFixed(decimals)}
              {suffix}
            </span>
          </div>
          <div className="stats-panel-stat-row">
            <span>Max:</span>
            <span>
              {stats.max.toFixed(decimals)}
              {suffix}
            </span>
          </div>
          <div className="stats-panel-stat-row">
            <span>Q1:</span>
            <span>
              {stats.q1.toFixed(decimals)}
              {suffix}
            </span>
          </div>
          <div className="stats-panel-stat-row">
            <span>Q3:</span>
            <span>
              {stats.q3.toFixed(decimals)}
              {suffix}
            </span>
          </div>
        </div>
      </div>
    );
  };

  const renderClusterStats = () => {
    if (!clusterStats) return null;

    // Calculate additional statistics from clusterData if available
    let mostFrequentSize = null;
    let avgClusterSize = null;

    if (clusterData && clusterData.length > 0) {
      // Count cluster sizes
      const clusterSizes: { [key: number]: number } = {};
      for (const item of clusterData) {
        clusterSizes[item.cluster] = (clusterSizes[item.cluster] || 0) + 1;
      }

      // Find most frequent cluster size
      const sizeCounts: { [key: number]: number } = {};
      for (const size of Object.values(clusterSizes)) {
        sizeCounts[size] = (sizeCounts[size] || 0) + 1;
      }

      let maxCount = 0;
      for (const [size, count] of Object.entries(sizeCounts)) {
        if (count > maxCount) {
          maxCount = count;
          mostFrequentSize = Number.parseInt(size);
        }
      }

      // Calculate average cluster size
      const totalItems = clusterData.length;
      const totalClusters = Object.keys(clusterSizes).length;
      avgClusterSize =
        totalClusters > 0 ? (totalItems / totalClusters).toFixed(1) : 0;
    }

    return (
      <div className="stats-panel-section">
        <div className="stats-panel-section-title">Cluster</div>
        <div className="stats-panel-section-content">
          <div className="stats-panel-stat-row">
            <span>Total Clusters:</span>
            <span>{clusterStats.total_clusters}</span>
          </div>
          <div className="stats-panel-stat-row">
            <span>Largest Cluster:</span>
            <span>{clusterStats.largest_cluster}</span>
          </div>
          <div className="stats-panel-stat-row">
            <span>Singletons:</span>
            <span>{clusterStats.singleton_clusters}</span>
          </div>
          {mostFrequentSize !== null && (
            <div className="stats-panel-stat-row">
              <span>Most Frequent Size:</span>
              <span>{mostFrequentSize}</span>
            </div>
          )}
          {avgClusterSize !== null && (
            <div className="stats-panel-stat-row">
              <span>Average Size:</span>
              <span>{avgClusterSize}</span>
            </div>
          )}
        </div>
      </div>
    );
  };

  const needsSeparator = () => {
    const hasAlignmentBar =
      isLzani && unalignedCount > 0 && activeDataSet === "scores";
    const hasStats = stats || clusterStats;
    return hasAlignmentBar && hasStats;
  };

  return (
    <>
      {/* Collapsed button */}
      {!isExpanded && (
        <button
          type="button"
          onClick={() => setIsExpanded(true)}
          className="stats-panel-button"
          title="Show statistics"
        >
          📊
        </button>
      )}

      {/* Expanded panel */}
      {isExpanded && (
        <div className="stats-panel">
          <div className="stats-panel-header">
            <h3>{getPanelTitle()}</h3>
            <button
              type="button"
              onClick={() => setIsExpanded(false)}
              title="Close"
            >
              ✕
            </button>
          </div>

          {/* Alignment bar for LZANI */}
          {renderAlignmentBar()}

          {/* Separator if needed */}
          {needsSeparator() && <div className="stats-panel-separator" />}

          {/* Statistics based on panel type */}
          {panelType !== "cluster" && renderDistributionStats()}
          {(panelType === "cluster" || panelType === "alignment") &&
            renderClusterStats()}
        </div>
      )}
    </>
  );
};
