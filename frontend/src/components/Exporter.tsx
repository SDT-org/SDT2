import React from "react";
import type { AppState, DocState, SetAppState } from "../appState";
import { useDocState } from "../hooks/useDocState";

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

  const runExport = React.useCallback(async () => {
    if (!tabs) return;

    try {
      // each tab exports its own images
      for (const tab of tabs) {
        currentTab.current = tab;
        if (tab !== docState.dataView) {
          swapDataView(tab);
        }
        await new Promise<void>((resolve) => {
          currentResolver.current = resolve;
        });
      }

      const result = await window.pywebview.api.export({
        doc_id: docState.id,
        export_path: appState.dataExportPath,
        output_cluster: true,
        cluster_threshold: docState.clustermap.threshold,
        cluster_method: docState.clustermap.method,
        image_format: appState.saveFormat,
        prefix: docState.exportPrefix,
        open_folder: docState.openExportFolder,
      });

      if (result) {
        setExportStatus("success");
      } else {
        setExportStatus("idle");
      }
      currentTab.current = null;
      swapDataView(initialTab.current);
    } catch (error) {
      setExportStatus("idle");
      currentTab.current = null;
      throw error;
    }
  }, [
    tabs,
    appState.dataExportPath,
    appState.saveFormat,
    docState.id,
    docState.clustermap.threshold,
    docState.clustermap.method,
    swapDataView,
    setExportStatus,
    docState.exportPrefix,
    docState.openExportFolder,
    docState.dataView,
  ]);

  React.useEffect(() => {
    if (appState.exportStatus === "exporting" && !currentTab.current) {
      initialTab.current = docState.dataView;
      runExport();
    }
  }, [appState.exportStatus, runExport, docState.dataView]);

  React.useEffect(() => {
    if (appState.exportStatus !== "success") {
      return;
    }

    const handler = () => {
      setExportStatus("idle");
      setAppState((prev) => ({
        ...prev,
        showExportModal: false,
      }));
    };

    const id = setTimeout(handler, 2000);

    return () => {
      clearTimeout(id);
    };
  }, [appState, setAppState, setExportStatus]);

  return (
    <RenderStatusContext.Provider value={{ markTabAsRendered, currentTab }}>
      {children}
    </RenderStatusContext.Provider>
  );
};
