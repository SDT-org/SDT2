import * as d3 from "d3";
import { useCallback, useEffect, useMemo, useState } from "react";
import { distinctColor } from "../colors";

export interface MetadataColorConfig {
  docId: string;
  columnName?: string | undefined;
  enabled: boolean;
}

export interface MetadataColorResult {
  values: Record<string, string | number> | null;
  getColor: (id: string) => string;
  columnType: "numeric" | "categorical" | null;
  loading: boolean;
  error: string | null;
  uniqueValues?: string[] | number[];
  colorScale?: d3.ScaleOrdinal<string, string> | d3.ScaleSequential<string>;
}

/**
 * Hook for handling metadata-based coloring across different visualizations
 * Designed to be reusable for UMAP, clustermap, heatmap, etc.
 */
export function useMetadataColors({
  docId,
  columnName,
  enabled,
}: MetadataColorConfig): MetadataColorResult {
  const [values, setValues] = useState<Record<string, string | number> | null>(
    null,
  );
  const [columnType, setColumnType] = useState<
    "numeric" | "categorical" | null
  >(null);
  const [loading, setLoading] = useState(false);
  const [error, setError] = useState<string | null>(null);

  // Fetch metadata when column changes
  useEffect(() => {
    if (!enabled || !columnName || !docId) {
      setValues(null);
      setColumnType(null);
      return;
    }

    const fetchMetadata = async () => {
      setLoading(true);
      setError(null);

      try {
        const response = await window.pywebview.api.data.get_metadata_for_umap(
          docId,
          columnName,
        );

        setValues(response.value_map);
        setColumnType(response.column_type);
      } catch (err) {
        console.error("Error fetching metadata:", err);
        setError(
          err instanceof Error ? err.message : "Failed to load metadata",
        );
        setValues(null);
        setColumnType(null);
      } finally {
        setLoading(false);
      }
    };

    fetchMetadata();
  }, [docId, columnName, enabled]);

  // Calculate unique values and color scale
  const { uniqueValues, colorScale } = useMemo(() => {
    if (!values || !columnType) {
      return { uniqueValues: undefined, colorScale: undefined };
    }

    const valueArray = Object.values(values);

    if (columnType === "categorical") {
      // Get unique string values
      const unique = Array.from(
        new Set(valueArray.filter((v): v is string => typeof v === "string")),
      );

      // Create ordinal scale with distinct colors
      const scale = d3
        .scaleOrdinal<string, string>()
        .domain(unique)
        .range(unique.map((_, i) => distinctColor(i + 1))); // +1 to avoid gray (cluster 0)

      return { uniqueValues: unique, colorScale: scale };
    }

    if (columnType === "numeric") {
      // Get numeric range
      const numericValues = valueArray.filter(
        (v): v is number => typeof v === "number",
      );
      const min = Math.min(...numericValues);
      const max = Math.max(...numericValues);

      // Create sequential scale (blue to red)
      const scale = d3
        .scaleSequential()
        .domain([min, max])
        .interpolator(d3.interpolateRdBu)
        .clamp(true);

      return {
        uniqueValues: [min, max], // Store range for legend
        colorScale: scale,
      };
    }

    return { uniqueValues: undefined, colorScale: undefined };
  }, [values, columnType]);

  // Color getter function
  const getColor = useCallback(
    (id: string): string => {
      if (!values || !colorScale) {
        return "#cccccc"; // Default gray
      }

      const value = values[id];
      if (value === undefined || value === null) {
        return "#cccccc"; // Gray for missing values
      }

      if (columnType === "categorical" && typeof value === "string") {
        return (colorScale as d3.ScaleOrdinal<string, string>)(value);
      }

      if (columnType === "numeric" && typeof value === "number") {
        return (colorScale as d3.ScaleSequential<string>)(value);
      }

      return "#cccccc"; // Fallback
    },
    [values, colorScale, columnType],
  );

  const result: MetadataColorResult = {
    values,
    getColor,
    columnType,
    loading,
    error,
  };

  if (uniqueValues !== undefined) {
    result.uniqueValues = uniqueValues;
  }

  if (colorScale !== undefined) {
    result.colorScale = colorScale;
  }

  return result;
}

/**
 * Helper function to generate tooltip content for metadata
 */
export function getMetadataTooltip(
  _id: string,
  columnName: string | undefined,
  value: string | number | undefined,
): string {
  if (!columnName || value === undefined || value === null) {
    return "";
  }
  return `<br><strong>${columnName}:</strong> ${value}`;
}
