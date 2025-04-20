import type React from "react";

export type TooltipStyle = "modern" | "minimal" | "glass";

interface TooltipData {
  x: number;
  y: number;
  value: number;
  xLabel: string;
  yLabel: string;
}

interface HeatmapTooltipProps {
  data: TooltipData;
  roundTo: number;
  style?: TooltipStyle;
}

export const tooltipStyles = {
  modern: {
    container: {
      background: "rgba(0, 0, 0, 0.8)",
      color: "white",
      padding: "8px 12px",
      borderRadius: "6px",
      fontSize: "14px",
      boxShadow: "0 4px 6px rgba(0, 0, 0, 0.1)",
      border: "none",
    },
    row: {
      marginBottom: "4px",
    },
    label: {
      color: "rgba(255, 255, 255, 0.7)",
      marginRight: "4px",
    },
  },
  minimal: {
    container: {
      background: "white",
      color: "black",
      padding: "8px",
      border: "1px solid #ddd",
      borderRadius: "4px",
      fontSize: "12px",
      boxShadow: "0 2px 4px rgba(0, 0, 0, 0.05)",
    },
    row: {
      marginBottom: "2px",
    },
    label: {
      color: "#666",
      marginRight: "4px",
    },
  },
  glass: {
    container: {
      background: "rgba(255, 255, 255, 0.9)",
      backdropFilter: "blur(8px)",
      padding: "10px 15px",
      borderRadius: "8px",
      border: "1px solid rgba(255, 255, 255, 0.3)",
      boxShadow: "0 4px 12px rgba(0, 0, 0, 0.1)",
      fontSize: "13px",
      color: "black",
    },
    row: {
      marginBottom: "6px",
    },
    label: {
      color: "#555",
      marginRight: "4px",
      fontWeight: 500,
    },
  },
} as const;

export const HeatmapTooltip: React.FC<HeatmapTooltipProps> = ({
  data,
  roundTo,
  style = "modern",
}) => {
  const styles = tooltipStyles[style];

  return (
    <div
      style={{
        position: "absolute",
        left: data.x + 10,
        top: data.y - 10,
        pointerEvents: "none",
        transform: "translate(0, -100%)",
        ...styles.container,
      }}
    >
      <div style={styles.row}>
        <span style={styles.label}>Row:</span>
        {data.xLabel}
      </div>
      <div style={styles.row}>
        <span style={styles.label}>Column:</span>
        {data.yLabel}
      </div>
      <div style={styles.row}>
        <span style={styles.label}>Value:</span>
        {data.value.toFixed(roundTo)}
      </div>
    </div>
  );
};
