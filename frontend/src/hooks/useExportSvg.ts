import React from "react";
import useAppState, { type SaveableImageKey } from "../appState";
import { useRenderStatus } from "../components/Exporter";
import { useDocState } from "./useDocState";

export const useExportSvg = (
  key: Extract<SaveableImageKey, "clustermap" | "heatmap">,
  ref: React.RefObject<SVGSVGElement>,
) => {
  const running = React.useRef(false);
  const { appState, setAppState } = useAppState();
  const { docState } = useDocState(
    appState.activeDocumentId,
    appState,
    setAppState,
  );
  const { markTabAsRendered, currentTab } = useRenderStatus();

  const exportSvg = React.useCallback(() => {
    if (
      currentTab.current !== key ||
      appState.exportStatus !== "exporting" ||
      !ref.current?.id
    ) {
      return;
    }

    const selector = `#${ref.current.id}`;

    const doExport = async () => {
      if (running.current) {
        return;
      }
      running.current = true;
      await window.pywebview.api.save_svg_element(docState.id, selector, key);
      markTabAsRendered(docState.dataView);
      running.current = false;
    };

    doExport();
  }, [
    appState,
    key,
    ref,
    currentTab,
    docState,
    appState.exportStatus,
    markTabAsRendered,
  ]);

  return exportSvg;
};
