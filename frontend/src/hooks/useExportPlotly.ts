import Plotly from "plotly.js-cartesian-dist-min";
import React from "react";
import useAppState, { type SaveableImageKey } from "../appState";
import { useRenderStatus } from "../components/Exporter";
import { assertDefined } from "../helpers";
import { useDocState } from "./useDocState";

// TODO: do this the good way if react-plotly gets ref support
const getPlotlyElement = () =>
  assertDefined(
    (
      document.getElementsByClassName(
        "js-plotly-plot",
      ) as HTMLCollectionOf<HTMLElement>
    )[0],
  );

export const useExportPlotly = (
  key: Extract<
    SaveableImageKey,
    "distribution_histogram" | "distribution_violin"
  >,
) => {
  const running = React.useRef(false);
  const { appState, setAppState } = useAppState();
  const { docState } = useDocState(
    appState.activeDocumentId,
    appState,
    setAppState,
  );

  const { markTabAsRendered, currentTab } = useRenderStatus();

  const exportPlotly = React.useCallback(() => {
    if (
      !(currentTab.current === key && appState.exportStatus === "exporting")
    ) {
      return;
    }

    const format = appState.saveFormat;

    const config = {
      format: appState.saveFormat,
      width: 1000,
      height: 800,
    };

    const doExport = async () => {
      if (running.current) {
        return;
      }
      running.current = true;

      await Plotly.toImage(getPlotlyElement(), config);
      const data = await Plotly.toImage(getPlotlyElement(), config);

      if (format === "svg") {
        await window.pywebview.api.export.save_svg_data(
          docState.id,
          data,
          key,
          format,
        );
      } else {
        await window.pywebview.api.export.save_raster_image(
          docState.id,
          data,
          key,
          format,
        );
      }

      markTabAsRendered(docState.dataView);
      running.current = false;
    };

    doExport();
  }, [
    key,
    appState.saveFormat,
    currentTab,
    docState,
    appState.exportStatus,
    markTabAsRendered,
  ]);

  return exportPlotly;
};
