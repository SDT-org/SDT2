import type { ClusterDataItem, HeatmapData } from "../plotTypes";
import type { DendrogramData } from "./dendrogram";

export interface ClustermapDataResponse {
  matrix: HeatmapData;
  tickText: string[];
  clusterData: ClusterDataItem[];
  dendrogramData?: DendrogramData | null;
}
