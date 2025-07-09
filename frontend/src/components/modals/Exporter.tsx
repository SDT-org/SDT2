import React from "react";
import { Button } from "react-aria-components";
import type { AppState, DocState, SetAppState } from "../../appState";
import { useDocState } from "../../hooks/useDocState";
import { SuccessModal } from "../ui/SuccessMessage";

type RenderStatus = {
  markTabAsRendered: (tabId: string) => void;
  currentTab: React.MutableRefObject<string | null>;
};

const RenderStatusContext = React.createContext<RenderStatus>({
  markTabAsRendered: () => {
    throw new Error("markTabAsRendered must be initialized");
  },
  currentTab: { current: null },
});

export const useRenderStatus = () => {
  const context = React.useContext(RenderStatusContext);
  if (!context) {
    throw new Error("useRenderStatus must be used within RenderStatusProvider");
  }
  return context;
};

interface ExporterProps {
  tabs: DocState["dataView"][];
  children: React.ReactNode;
  appState: AppState;
  setAppState: SetAppState;
  docState: DocState | undefined;
}

export const Exporter = ({
  tabs,
  children,
  appState,
  setAppState,
  docState,
}: ExporterProps) => {
  if (!docState) {
    return null;
  }
  const { updateDocState } = useDocState(docState.id, appState, setAppState);
  const initialTab = React.useRef<DocState["dataView"]>(docState.dataView);
  const currentTab = React.useRef<string | null>(null);
  const currentResolver = React.useRef<(() => void) | null>(null);
  const exportStep = React.useRef("Preparing");

  const swapDataView = React.useCallback(
    (view: DocState["dataView"]) => {
      updateDocState({ dataView: view });
    },
    [updateDocState],
  );

  const setExportStatus = React.useCallback(
    (status: AppState["exportStatus"]) =>
      setAppState((prev) => ({ ...prev, exportStatus: status })),
    [setAppState],
  );

  const markTabAsRendered = React.useCallback((tabId: string) => {
    if (!currentTab.current) return;

    if (tabId === currentTab.current && currentResolver.current) {
      currentResolver.current();
    } else {
      console.warn(
        `tried completing ${tabId}, but currently processing ${currentTab.current}`,
      );
    }
  }, []);

  const cancelExport = React.useCallback(() => {
    if (currentTab.current) {
      currentTab.current = null;
      currentResolver.current = null;
      exportStep.current = "Preparing";
      setExportStatus("idle");
      updateDocState({ dataView: initialTab.current });
    }
  }, [setExportStatus, updateDocState]);

  const setFinalStage = React.useCallback(
    (stage: AppState["exportStatus"]) => {
      setExportStatus(stage);
      exportStep.current = "Preparing";
      currentTab.current = null;
      currentResolver.current = null;
    },
    [setExportStatus],
  );

  const runExport = React.useCallback(async () => {
    if (currentTab.current || !tabs) return;

    setExportStatus("exporting");

    try {
      // each tab exports its own images
      for (const tab of tabs) {
        exportStep.current = `Exporting ${tab}`;
        currentTab.current = tab;

        swapDataView(tab);

        await new Promise<void>((resolve) => {
          currentResolver.current = resolve;
        });
      }

      currentTab.current = null;
      currentResolver.current = null;
      exportStep.current = "Finalizing";
      swapDataView(initialTab.current);

      const result = await window.pywebview.api.export({
        doc_id: docState.id,
        export_path: appState.dataExportPath,
        output_cluster: true,
        cluster_threshold: docState.clustermap.threshold,
        cluster_method: docState.clustermap.method,
        image_format: appState.saveFormat,
        prefix: docState.exportPrefix,
        open_folder: appState.openExportFolder,
      });

      if (result) {
        setFinalStage("success");
      } else {
        setFinalStage("idle");
      }
    } catch (error) {
      setFinalStage("idle");
      setAppState((prev) => ({
        ...prev,
        error: error as Error,
      }));
    }
  }, [
    tabs,
    appState.dataExportPath,
    appState.saveFormat,
    appState.openExportFolder,
    docState.id,
    docState.clustermap.threshold,
    docState.clustermap.method,
    swapDataView,
    setExportStatus,
    docState.exportPrefix,
    setFinalStage,
    setAppState,
  ]);

  React.useEffect(() => {
    if (
      appState.exportStatus === "preparing" &&
      !currentTab.current &&
      currentResolver.current === null
    ) {
      runExport();
    }
  }, [runExport, appState.exportStatus]);

  React.useEffect(() => {
    if (appState.exportStatus !== "success") {
      return;
    }

    const handler = () => {
      setExportStatus("idle");
    };

    const id = setTimeout(handler, 3000);

    return () => {
      clearTimeout(id);
    };
  }, [appState, setExportStatus]);

  return (
    <RenderStatusContext.Provider value={{ markTabAsRendered, currentTab }}>
      {appState.exportStatus !== "idle" ? (
        appState.exportStatus === "success" ? (
          <SuccessModal onClose={() => setExportStatus("idle")} />
        ) : (
          <div className="app-backdrop">
            <div className="app-overlay">
              <div className="app-loader" />
              <div className="app-loader-status">
                {currentTab.current
                  ? `Exporting ${currentTab.current?.replaceAll("_", " ")}`
                  : exportStep.current}
                ...
              </div>
            </div>
            <Button
              className="react-aria-Button app-overlay-cancel"
              onPress={cancelExport}
            >
              Cancel
            </Button>
          </div>
        )
      ) : null}
      {children}
    </RenderStatusContext.Provider>
  );
};
