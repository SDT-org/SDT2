import React from "react";
import { type DocState, type SetDocState, initialDocState } from "../appState";
import type {
  DistributionData,
  GetDataResponse,
  HeatmapData,
} from "../plotTypes";
import { getScaledFontSize } from "../plotUtils";
import { services } from "../services";

export const useGetData = (docState: DocState, setDocState: SetDocState) => {
  const [loading, setLoading] = React.useState(false);
  const [tickText, setTickText] = React.useState<string[]>([""]);
  const [heatmapData, setHeatmapData] = React.useState<HeatmapData>();
  const [distributionData, setDistributionData] =
    React.useState<DistributionData>();
  const [metaData, setMetaData] = React.useState<
    GetDataResponse["metadata"] | undefined
  >();

  React.useEffect(() => {
    setLoading(true);

    services
      .getData(docState.id)
      .then((rawData) => {
        const parsedResponse: GetDataResponse = JSON.parse(
          rawData.replace(/\bNaN\b/g, "null"),
        );
        const { data, metadata, identity_scores, gc_stats, length_stats } =
          parsedResponse;

        const [tickText, ...parsedData] = data;

        setMetaData(metadata);
        setTickText(tickText as string[]);
        setHeatmapData(parsedData);
        setDistributionData({
          gc_stats,
          length_stats,
          raw_mat: identity_scores.map((i) => i[2]),
          identity_combos: identity_scores.map((i) => [i[0], i[1]]),
        });

        setDocState((prev) => ({
          ...prev,
          heatmap: {
            ...prev.heatmap,
            vmin: metadata.minVal,
            ...(docState.sequences_count > 99 && {
              annotation: false,
              axis_labels: false,
              cellspace: 0,
            }),
            annotation_font_size: getScaledFontSize(
              initialDocState.heatmap.annotation_font_size,
              parsedData.length,
            ),
          },
        }));
      })
      .finally(() => {
        setLoading(false);
      });
  }, [docState.id, docState.sequences_count, setDocState]);

  return {
    loading,
    tickText,
    heatmapData,
    distributionData,
    metaData,
  };
};
