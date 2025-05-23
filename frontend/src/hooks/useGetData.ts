import React from "react";
import { type DocState, type SetDocState, initialDocState } from "../appState";
import type {
  DistributionData,
  GetDataResponse,
  HeatmapData,
} from "../plotTypes";
import { getScaledFontSize } from "../plotUtils";
import { services } from "../services";
import { getDocument } from "../services/documents";

export const useGetData = (docState: DocState, setDocState: SetDocState) => {
  const [loading, setLoading] = React.useState(false);
  const [tickText, setTickText] = React.useState<string[]>([""]);
  const [heatmapData, setHeatmapData] = React.useState<HeatmapData>();
  const [distributionData, setDistributionData] =
    React.useState<DistributionData>();
  const [metaData, setMetaData] = React.useState<
    GetDataResponse["metadata"] | undefined
  >();

  // TODO: replace this entire thing, docstate already gets pushed by the backend
  // so only make this set the defaults. also move all the data into the document
  // so we aren't parsing it out here.
  React.useEffect(() => {
    setLoading(true);
    services
      .getData(docState.id)
      .then(async (rawData) => {
        const parsedResponse: GetDataResponse = JSON.parse(
          rawData.replace(/\bNaN\b/g, "null"),
        );

        const { data, metadata, ids, identity_scores, full_stats } =
          parsedResponse;
        const [tickText, ...parsedData] = data;

        setMetaData(metadata);
        setTickText(tickText as string[]);
        setHeatmapData(parsedData);
        setDistributionData({
          full_stats,
          gc_stats: full_stats.map((row) => row[1]),
          length_stats: full_stats.map((row) => row[2]),
          raw_mat: identity_scores.map((i) => i[2]),
          ids,
          identity_combos: identity_scores.map((i) => [i[0], i[1]]),
        });

        const state = await getDocument(docState.id);
        const scaledAxisLabelFontSize = getScaledFontSize(
          initialDocState.heatmap.axlabel_fontsize,
          parsedData.map(Boolean).length,
        );

        setDocState(
          (prev) => ({
            ...prev,
            ...state,
            heatmap: {
              ...prev.heatmap,
              ...state.heatmap,
              vmin: metadata.minVal,
              ...(docState.sequences_count > 99
                ? {
                    annotation: false,
                    axis_labels: false,
                    cellspace: 0,
                  }
                : null),
              ...(docState.filetype === "application/vnd.sdt"
                ? null
                : {
                    axlabel_fontsize: scaledAxisLabelFontSize,
                    axlabel_yfontsize: scaledAxisLabelFontSize,
                  }),
            },
            clustermap: {
              ...prev.clustermap,
              ...state.clustermap,
              ...(docState.filetype === "application/vnd.sdt"
                ? null
                : {
                    axlabel_fontsize: scaledAxisLabelFontSize,
                    axlabel_yfontsize: scaledAxisLabelFontSize,
                  }),
            },
            ...(docState.sequences_count > 99
              ? {
                  distribution: {
                    ...prev.distribution,
                    ...state.distribution,
                    violin: {
                      ...prev.distribution.violin,
                      ...state.distribution.violin,
                      showPoints: false,
                    },
                  },
                }
              : null),
          }),
          false,
        );
      })
      .finally(() => {
        setLoading(false);
      });
  }, [docState.id, docState.filetype, docState.sequences_count, setDocState]);

  return {
    loading,
    tickText,
    heatmapData,
    distributionData,
    metaData,
  };
};
