import React from "react";
import { type Key, Tab, TabList, TabPanel, Tabs } from "react-aria-components";
import { type AppState, AppStateContext, initialAppState } from "../appState";
import { useAppBlur } from "../hooks/appBlur";
import { useWaitForPywebview } from "../hooks/usePywebviewReadyEvent";
import { useSaveState } from "../hooks/useSaveState";
import { useShortcutKeys } from "../hooks/useShortcutKeys";
import { useSyncState } from "../hooks/useSyncState";
import { Document } from "./Document";
// import { restoreClientState } from "../restoreClientState";
import { ErrorBoundary } from "./ErrorBoundary";
import { ExportModal } from "./ExportModal";

export const App = () => {
  // const restoreInitialState = React.useCallback(() => {
  //   setLoading(true);

  //   window.pywebview.api.get_state().then((data) => {
  //     setAppState((prev) => ({
  //       ...data,
  //       client: restoreClientState(prev.client),
  //       showExportModal: false,
  //     }));
  //     setLoading(false);
  //     setInitialized(true);
  //   });
  // }, []);

  const [appState, setAppState] = React.useState<AppState>(initialAppState);
  // const [loading, setLoading] = React.useState(true);
  // const [loading] = React.useState(true);
  const [initialized, setInitialized] = React.useState(false);
  // const startRun = useStartRun(appState);
  useSyncState(setAppState);
  useShortcutKeys(appState, setAppState);
  useAppBlur();
  useWaitForPywebview(() => setInitialized(true));
  // useWaitForPywebview(restoreInitialState);
  useSaveState(initialized, appState);

  React.useEffect(() => {
    if (!initialized) {
      return;
    }
    if (appState.documents.length === 0) {
      window.pywebview.api.new_doc();
      return;
    }

    if (
      appState.documents.length > 0 &&
      !appState.documents.find((d) => d.id === appState.activeDocumentId)
    ) {
      setAppState((prev) => {
        const newActive = appState.documents[appState.documents.length - 1];
        if (!newActive) {
          return prev;
        }
        return {
          ...prev,
          activeDocumentId: newActive.id,
        };
      });
    }
  }, [initialized, appState.documents, appState.activeDocumentId]);

  const setActiveDocumentId = (id: Key) =>
    setAppState((prev) => ({ ...prev, activeDocumentId: id as string }));

  const tabView = "tabs";

  return (
    <ErrorBoundary appState={appState} setAppState={setAppState}>
      <AppStateContext.Provider value={{ appState, setAppState }}>
        {initialized &&
        appState.documents.length > 0 &&
        appState.activeDocumentId ? (
          <>
            <Tabs
              selectedKey={appState.activeDocumentId}
              onSelectionChange={setActiveDocumentId}
            >
              {tabView === "tabs" ? (
                <TabList>
                  {appState.documents.map((doc) => (
                    <Tab id={doc.id} key={doc.id}>
                      {doc.basename || "Untitled"}
                      {doc.stage === "Analyzing" ? (
                        ` (${doc.progress}%)`
                      ) : doc.stage === "Postprocessing" ? (
                        <div className="app-loader" />
                      ) : null}
                    </Tab>
                  ))}
                </TabList>
              ) : null}
              {appState.documents.map((doc) => (
                <TabPanel
                  id={doc.id}
                  key={doc.id}
                  shouldForceMount={doc.view === "loader"}
                  style={{
                    display:
                      doc.id === appState.activeDocumentId ? "block" : "none",
                  }}
                >
                  <Document id={doc.id} key={doc.id} tabView={tabView} />
                </TabPanel>
              ))}
            </Tabs>
            {appState.activeDocumentId ? <ExportModal /> : null}
          </>
        ) : (
          <div className="app-overlay app-loader" />
        )}
      </AppStateContext.Provider>
    </ErrorBoundary>
  );
};
