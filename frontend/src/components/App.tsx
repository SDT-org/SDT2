// import Plotly from "plotly.js-cartesian-dist-min";
import React from "react";
import { Button, type Key, Tab, TabList, Tabs } from "react-aria-components";
import {
  TbFile,
  TbFolderShare,
  TbLayoutSidebarLeftCollapse,
  TbLayoutSidebarLeftExpand,
  TbPlus,
} from "react-icons/tb";
import {
  type AppState,
  AppStateContext,
  findDoc,
  initialAppState,
} from "../appState";
import { useAppBlur } from "../hooks/appBlur";
import { useCloseDocument } from "../hooks/useCloseDocument";
import useNewDocument from "../hooks/useNewDocument";
import { useWaitForPywebview } from "../hooks/usePywebviewReadyEvent";
import { useShortcutKeys } from "../hooks/useShortcutKeys";
import { useSyncState } from "../hooks/useSyncState";
import { Document } from "./Document";
import { ErrorBoundary } from "./ErrorBoundary";
import { ExportModal } from "./ExportModal";
import { Exporter } from "./Exporter";
import { MainMenu } from "./Menu";
import { Select, SelectItem } from "./Select";

export const App = () => {
  const [appState, setAppState] = React.useState<AppState>(initialAppState);
  const [initialized, setInitialized] = React.useState(false);

  const handleReady = React.useCallback(() => {
    if (initialized) {
      return;
    }

    Promise.all([
      window.pywebview.api.system.get_app_config(),
      window.pywebview.api.system.get_app_settings(),
    ]).then(
      ([
        config,
        {
          user_settings: { export_path, open_folder_after_export },
        },
      ]) => {
        setAppState((prev) => ({
          ...prev,
          config,
          dataExportPath: export_path,
          openExportFolder:
            open_folder_after_export || initialAppState.openExportFolder,
        }));
        setInitialized(true);
      },
    );
  }, [initialized]);

  useSyncState(setAppState);
  useShortcutKeys(appState, setAppState);
  useAppBlur();
  useWaitForPywebview(handleReady);
  const tabListRef = React.useRef<HTMLDivElement>(null);
  const lastTabRef = React.useRef<HTMLDivElement>(null);

  React.useEffect(() => {
    if (!initialized) {
      return;
    }
    if (appState.documents.length === 0) {
      window.pywebview.api.documents.new_doc();
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

  React.useEffect(() => {
    const id = setTimeout(() => {
      if (!tabListRef.current || !appState.activeDocumentId) {
        return;
      }
      const activeTab = tabListRef.current.querySelector("[data-selected]");

      if (activeTab) {
        activeTab?.scrollIntoView();
      }

      const doc = findDoc(appState.activeDocumentId, appState.documents);
      if (doc) {
        window.pywebview.api.documents.set_window_title(
          doc.basename || "Untitled",
        );
      }
    }, 0);

    return () => clearTimeout(id);
  }, [appState.activeDocumentId, appState.documents]);

  const setActiveDocumentId = React.useCallback((id: Key | null) => {
    if (id !== null) {
      setAppState((prev) => ({ ...prev, activeDocumentId: id as string }));
    }
  }, []);

  const tabView = React.useMemo(
    () => (appState.documents.length > 20 ? "select" : "tabs"),
    [appState.documents.length],
  );

  const closeDocument = useCloseDocument(appState, setAppState);
  const newDocument = useNewDocument(setAppState);

  const [leftSidebarCollapsed, setLeftSidebarCollapsed] = React.useState(false);
  const [newTabPositionRect, setNewTabPositionRect] = React.useState<DOMRect>();
  const [hideNewTab, setHideNewTab] = React.useState(true);

  const checkTabOverflow = React.useCallback(() => {
    const container = tabListRef.current;
    if (!container) {
      return false;
    }

    return (
      container.scrollWidth > container.offsetWidth ||
      container.scrollHeight > container.offsetHeight
    );
  }, []);

  const updateNewTabOffset = React.useCallback(() => {
    if (!lastTabRef.current) {
      return;
    }
    const rect = checkTabOverflow()
      ? tabListRef?.current?.getBoundingClientRect()
      : lastTabRef.current.getBoundingClientRect();

    setNewTabPositionRect(rect);
    setHideNewTab(false);
  }, [checkTabOverflow]);

  React.useEffect(() => {
    if (!appState.documents.length) {
      return;
    }
    const id = setTimeout(() => updateNewTabOffset(), 0);
    return () => clearTimeout(id);
  }, [updateNewTabOffset, appState.documents.length]);

  React.useEffect(() => {
    let id: Timer;

    const observer = new MutationObserver(() => {
      id = setTimeout(updateNewTabOffset, 0);
    });

    if (tabListRef.current) {
      observer.observe(tabListRef.current, {
        attributes: false,
        childList: true,
        subtree: true,
      });
    }

    window.addEventListener("resize", updateNewTabOffset);

    return () => {
      clearTimeout(id);
      window.removeEventListener("resize", updateNewTabOffset);
      observer.disconnect();
    };
  }, [updateNewTabOffset]);

  const activeDocument = React.useMemo(
    () =>
      appState.documents.find((doc) => doc.id === appState.activeDocumentId),
    [appState.activeDocumentId, appState.documents],
  );

  return (
    <ErrorBoundary appState={appState} setAppState={setAppState}>
      <AppStateContext.Provider value={{ appState, setAppState }}>
        {initialized &&
        appState.documents.length > 0 &&
        appState.activeDocumentId ? (
          <Exporter
            appState={appState}
            setAppState={setAppState}
            docState={activeDocument}
            tabs={[
              "heatmap",
              "clustermap",
              "distribution_histogram",
              "distribution_violin",
            ]}
          >
            <div className="app-layout">
              <div className="app-header">
                <div className="app-header-left">
                  <MainMenu />
                  {activeDocument && activeDocument.view === "viewer" ? (
                    <Button
                      onPress={() => setLeftSidebarCollapsed((prev) => !prev)}
                    >
                      <span data-sr-only>
                        {leftSidebarCollapsed
                          ? "Expand sidebar"
                          : "Collapse sidebar"}
                      </span>
                      {leftSidebarCollapsed ? (
                        <TbLayoutSidebarLeftExpand size={18} />
                      ) : (
                        <TbLayoutSidebarLeftCollapse size={18} />
                      )}
                    </Button>
                  ) : null}
                </div>
                <div className="app-header-middle">
                  {tabView === "tabs" ? (
                    <Tabs
                      selectedKey={appState.activeDocumentId}
                      onSelectionChange={setActiveDocumentId}
                      className="document-tabs"
                    >
                      <TabList className="document-tablist" ref={tabListRef}>
                        {appState.documents.map((doc, index) => (
                          <Tab
                            id={doc.id}
                            key={doc.id}
                            className="react-aria-Tab document-tab"
                            ref={
                              appState.documents.length - 1 === index
                                ? lastTabRef
                                : null
                            }
                          >
                            <span
                              className="tab-title"
                              data-modified={doc.modified}
                            >
                              {doc.basename || "Untitled"}
                              <i>{doc.modified ? "Edited" : null}</i>
                            </span>
                            {doc.stage === "Analyzing" ? (
                              <span className="tab-progress text-tabular-nums">
                                {" "}
                                {doc.progress}%
                              </span>
                            ) : (
                              <>
                                {!doc.stage ? (
                                  <Button
                                    className={"close-button"}
                                    aria-label={`Close ${doc.basename}`}
                                    onPress={() => closeDocument(doc.id)}
                                  >
                                    <svg
                                      aria-hidden={true}
                                      xmlns="http://www.w3.org/2000/svg"
                                      viewBox="0 0 24 24"
                                      width="9"
                                      height="9"
                                      fill="none"
                                      stroke="currentColor"
                                      strokeWidth="4"
                                      strokeLinecap="round"
                                      strokeLinejoin="round"
                                    >
                                      <line x1="4" y1="4" x2="20" y2="20" />
                                      <line x1="20" y1="4" x2="4" y2="20" />
                                    </svg>
                                  </Button>
                                ) : null}
                              </>
                            )}
                            <TbFile size={16} />
                          </Tab>
                        ))}
                      </TabList>
                      <Button
                        aria-label="New Document"
                        className={"react-aria-Button new-document"}
                        onPress={newDocument}
                        style={{
                          left: newTabPositionRect?.right || 0,
                          top: newTabPositionRect?.top,
                          opacity: hideNewTab ? "0" : "1",
                        }}
                      >
                        <TbPlus size={16} />
                      </Button>
                    </Tabs>
                  ) : (
                    <Select
                      wide
                      selectedKey={appState.activeDocumentId}
                      onSelectionChange={(value) =>
                        value && setActiveDocumentId(value)
                      }
                      items={appState.documents.map((doc) => ({
                        id: doc.id,
                        name: doc.basename,
                      }))}
                    >
                      {({ id, name }) => (
                        <SelectItem id={id} textValue={name}>
                          {name || "Untitled"}
                        </SelectItem>
                      )}
                    </Select>
                  )}
                </div>
                {activeDocument?.view === "viewer" &&
                !activeDocument.invalid ? (
                  <div className="app-header-right">
                    <Button
                      className={"react-aria-Button with-svg"}
                      onPress={() =>
                        setAppState((prev) => ({
                          ...prev,
                          showExportModal: true,
                        }))
                      }
                    >
                      <TbFolderShare size={17} />
                      <span>Export</span>
                    </Button>
                  </div>
                ) : null}
              </div>
              <Tabs
                selectedKey={appState.activeDocumentId}
                onSelectionChange={setActiveDocumentId}
                className={"app-body"}
              >
                {appState.documents
                  .filter((doc) => doc.parsed)
                  .map((doc) => (
                    <Document
                      id={doc.id}
                      key={doc.id}
                      tabView={tabView}
                      leftSidebarCollapsed={leftSidebarCollapsed}
                    />
                  ))}
              </Tabs>
              {activeDocument ? (
                <ExportModal appState={appState} setAppState={setAppState} />
              ) : null}
            </div>
          </Exporter>
        ) : (
          <div className="app-overlay app-loader" />
        )}
      </AppStateContext.Provider>
    </ErrorBoundary>
  );
};
