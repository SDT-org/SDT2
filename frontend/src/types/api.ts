import type { ClusterDataItem, HeatmapData } from "../plotTypes";

export interface ClustermapDataResponse {
  matrix: HeatmapData;
  tickText: string[];
  clusterData: ClusterDataItem[];
}
