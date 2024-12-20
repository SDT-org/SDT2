import React from "react";
import useAppState, { initialAppState } from "../appState";
import type {
  DistributionData,
  GetDataResponse,
  HeatmapData,
} from "../plotTypes";
import { getScaledFontSize } from "../plotUtils";
import { services } from "../services";

export const useGetData = () => {
  const { appState, setAppState } = useAppState();
  const [loading, setLoading] = React.useState(false);
  const [tickText, setTickText] = React.useState<string[]>([""]);
  const [heatmapData, setHeatmapData] = React.useState<HeatmapData>();
  const [distributionData, setDistributionData] =
    React.useState<DistributionData>();
  const [metaData, setMetaData] = React.useState<
    GetDataResponse["metadata"] | undefined
  >();

  const getData = React.useCallback(() => {
    setLoading(true);

    services
      .getData()
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
        setAppState((prev) => ({
          ...prev,
          client: {
            ...prev.client,
            heatmap: {
              ...prev.client.heatmap,
              vmin: metadata.minVal,
              ...(appState.sequences_count > 99 && {
                annotation: false,
                axis_labels: false,
                cellspace: 0,
              }),
              annotation_font_size: getScaledFontSize(
                initialAppState.client.heatmap.annotation_font_size,
                parsedData.length,
              ),
            },
          },
        }));
      })
      .finally(() => {
        setLoading(false);
      });
  }, [appState.sequences_count, setAppState]);

  React.useEffect(() => {
    if (!appState.filemtime) {
      return;
    }
    getData();
  }, [getData, appState.filemtime]);

  return {
    getData,
    loading,
    tickText,
    heatmapData,
    distributionData,
    metaData,
  };
};
