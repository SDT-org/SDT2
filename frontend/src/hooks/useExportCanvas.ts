import React from "react";
import useAppState, { type SaveableImageKey } from "../appState";
import { useRenderStatus } from "../components/modals/Exporter";
import { isRasterFormat } from "../helpers";
import { useDocState } from "./useDocState";

export const useExportCanvas = (
  key: Extract<SaveableImageKey, "clustermap" | "heatmap">,
  ref: React.RefObject<HTMLCanvasElement>,
) => {
  const running = React.useRef(false);
  const { appState, setAppState } = useAppState();
  const { docState } = useDocState(
    appState.activeDocumentId,
    appState,
    setAppState,
  );

  const { markTabAsRendered, currentTab } = useRenderStatus();

  const exportRaster = React.useCallback(() => {
    if (
      !isRasterFormat(appState.saveFormat) ||
      currentTab.current !== key ||
      appState.exportStatus !== "exporting" ||
      !ref.current
    ) {
      return;
    }

    const format = appState.saveFormat;

    const doExport = async () => {
      if (running.current) {
        return;
      }
      running.current = true;
      const data: string = await new Promise((resolve, reject) => {
        if (!ref.current) {
          throw new Error("Expected ref have a current value");
        }

        // work around slow canvas issues
        let attempts = 0;
        const maxAttempts = 5;
        const interval = setInterval(() => {
          if (attempts >= maxAttempts) {
            clearInterval(interval);
            reject(new Error("Blob conversion failed after 5 attempts"));
          } else {
            (ref.current as HTMLCanvasElement).toBlob(async (blob) => {
              if (blob) {
                clearInterval(interval);
                const arrayBuffer = await blob.arrayBuffer();
                const binary = Array.from(new Uint8Array(arrayBuffer))
                  .map((byte) => String.fromCharCode(byte))
                  .join("");
                resolve(`data:image/${format};base64,${btoa(binary)}`);
              }
            }, `image/${format}`);
            attempts++;
          }
        }, 60);
      });

      await window.pywebview.api.save_raster_image(
        docState.id,
        data,
        key,
        format,
      );

      markTabAsRendered(docState.dataView);
      running.current = false;
    };

    doExport();
  }, [
    key,
    ref,
    appState.saveFormat,
    currentTab,
    docState,
    appState.exportStatus,
    markTabAsRendered,
  ]);

  return exportRaster;
};
