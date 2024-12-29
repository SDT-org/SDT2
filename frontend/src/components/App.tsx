import React from "react";
import { Button, type Key, Tab, TabList, Tabs } from "react-aria-components";
import { type AppState, AppStateContext, initialAppState } from "../appState";
import { useAppBlur } from "../hooks/appBlur";
import { useCloseDocument } from "../hooks/useCloseDocument";
import { useWaitForPywebview } from "../hooks/usePywebviewReadyEvent";
import { useSaveState } from "../hooks/useSaveState";
import { useShortcutKeys } from "../hooks/useShortcutKeys";
import { useSyncState } from "../hooks/useSyncState";
import { Document } from "./Document";
// import { restoreClientState } from "../restoreClientState";
import { ErrorBoundary } from "./ErrorBoundary";
import { ExportModal } from "./ExportModal";
import { MainMenu } from "./Menu";

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
  // useSaveState(initialized, appState);

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
  const closeDocument = useCloseDocument(appState, setAppState);

  const [leftSidebarCollapsed, setLeftSidebarCollapsed] = React.useState(false);

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
              className={"react-aria-Tabs app-layout"}
            >
              {tabView === "tabs" ? (
                <div className="document-tabs appearance-alt appearance-alt-light appearance-alt-round">
                  <div className="sidebar-left-buttons">
                    <MainMenu />
                    <Button
                      onPress={() => setLeftSidebarCollapsed((prev) => !prev)}
                    >
                      <span data-sr-only>
                        {leftSidebarCollapsed
                          ? "Expand sidebar"
                          : "Collapse sidebar"}
                      </span>
                      <svg
                        xmlns="http://www.w3.org/2000/svg"
                        height="1.6rem"
                        fill="currentColor"
                        viewBox="0 0 512 512"
                        aria-hidden="true"
                      >
                        <path d="M 448 64 Q 462 64 471 73 L 471 73 L 471 73 Q 480 82 480 96 L 480 416 L 480 416 Q 480 430 471 439 Q 462 448 448 448 L 224 448 L 224 448 L 224 64 L 224 64 L 448 64 L 448 64 Z M 64 64 L 192 64 L 64 64 L 192 64 L 192 448 L 192 448 L 64 448 L 64 448 Q 50 448 41 439 Q 32 430 32 416 L 32 96 L 32 96 Q 32 82 41 73 Q 50 64 64 64 L 64 64 Z M 64 32 Q 37 33 19 51 L 19 51 L 19 51 Q 1 69 0 96 L 0 416 L 0 416 Q 1 443 19 461 Q 37 479 64 480 L 448 480 L 448 480 Q 475 479 493 461 Q 511 443 512 416 L 512 96 L 512 96 Q 511 69 493 51 Q 475 33 448 32 L 64 32 L 64 32 Z M 80 96 Q 65 97 64 112 Q 65 127 80 128 L 144 128 L 144 128 Q 159 127 160 112 Q 159 97 144 96 L 80 96 L 80 96 Z M 64 176 Q 65 191 80 192 L 144 192 L 144 192 Q 159 191 160 176 Q 159 161 144 160 L 80 160 L 80 160 Q 65 161 64 176 L 64 176 Z M 80 224 Q 65 225 64 240 Q 65 255 80 256 L 144 256 L 144 256 Q 159 255 160 240 Q 159 225 144 224 L 80 224 L 80 224 Z" />
                      </svg>
                    </Button>
                  </div>
                  {appState.documents.length > 0 ? (
                    <TabList className="react-aria-TabList document-tablist">
                      {appState.documents.map((doc) => (
                        <Tab id={doc.id} key={doc.id}>
                          <span>{doc.basename || "Untitled"}</span>
                          {doc.stage === "Analyzing" ? (
                            <span className="text-tabular-nums">
                              {" "}
                              {doc.progress}%
                            </span>
                          ) : doc.stage === "Postprocessing" ? (
                            <div className="app-loader" />
                          ) : (
                            <Button
                              className={"close-button"}
                              aria-label={`Close ${doc.basename}`}
                              onPress={() => closeDocument(doc.id)}
                            >
                              <svg
                                aria-hidden={true}
                                xmlns="http://www.w3.org/2000/svg"
                                viewBox="0 0 24 24"
                                width="10"
                                height="10"
                                fill="none"
                                stroke="currentColor"
                                strokeWidth="2"
                                strokeLinecap="round"
                                strokeLinejoin="round"
                              >
                                <line x1="4" y1="4" x2="20" y2="20" />
                                <line x1="20" y1="4" x2="4" y2="20" />
                              </svg>
                            </Button>
                          )}
                        </Tab>
                      ))}
                    </TabList>
                  ) : null}
                </div>
              ) : null}
              {appState.documents.map((doc) => (
                <Document
                  id={doc.id}
                  key={doc.id}
                  tabView={tabView}
                  leftSidebarCollapsed={leftSidebarCollapsed}
                />
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
