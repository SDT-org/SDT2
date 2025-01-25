import React from "react";
import { Button, type Key, Tab, TabList, Tabs } from "react-aria-components";
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
import { HeatmapRefProvider } from "./HeatmapRefProvider";
import Icons from "./Icons";
import { MainMenu } from "./Menu";
import { Select, SelectItem } from "./Select";

export const App = () => {
  const [appState, setAppState] = React.useState<AppState>(initialAppState);
  const [initialized, setInitialized] = React.useState(false);
  useSyncState(setAppState);
  useShortcutKeys(appState, setAppState);
  useAppBlur();
  useWaitForPywebview(() => {
    if (initialized) {
      return;
    }
    window.pywebview.api.app_config().then((data) => {
      setAppState((prev) => ({ ...prev, config: data }));
      setInitialized(true);
    });
  });
  const tabListRef = React.useRef<HTMLDivElement>(null);
  const lastTabRef = React.useRef<HTMLDivElement>(null);

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
        window.pywebview.api.set_window_title(doc.basename || "Untitled");
      }
    }, 0);

    return () => clearTimeout(id);
  }, [appState.activeDocumentId, appState.documents]);

  const setActiveDocumentId = React.useCallback(
    (id: Key) =>
      setAppState((prev) => ({ ...prev, activeDocumentId: id as string })),
    [],
  );

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
        <HeatmapRefProvider>
          {initialized &&
          appState.documents.length > 0 &&
          appState.activeDocumentId ? (
            <div className="app-layout">
              <div className="app-header">
                <div className="app-header-left">
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
                      height="16"
                      fill="currentColor"
                      viewBox="0 0 512 512"
                      aria-hidden="true"
                    >
                      <path d="M 448 64 Q 462 64 471 73 L 471 73 L 471 73 Q 480 82 480 96 L 480 416 L 480 416 Q 480 430 471 439 Q 462 448 448 448 L 224 448 L 224 448 L 224 64 L 224 64 L 448 64 L 448 64 Z M 64 64 L 192 64 L 64 64 L 192 64 L 192 448 L 192 448 L 64 448 L 64 448 Q 50 448 41 439 Q 32 430 32 416 L 32 96 L 32 96 Q 32 82 41 73 Q 50 64 64 64 L 64 64 Z M 64 32 Q 37 33 19 51 L 19 51 L 19 51 Q 1 69 0 96 L 0 416 L 0 416 Q 1 443 19 461 Q 37 479 64 480 L 448 480 L 448 480 Q 475 479 493 461 Q 511 443 512 416 L 512 96 L 512 96 Q 511 69 493 51 Q 475 33 448 32 L 64 32 L 64 32 Z M 80 96 Q 65 97 64 112 Q 65 127 80 128 L 144 128 L 144 128 Q 159 127 160 112 Q 159 97 144 96 L 80 96 L 80 96 Z M 64 176 Q 65 191 80 192 L 144 192 L 144 192 Q 159 191 160 176 Q 159 161 144 160 L 80 160 L 80 160 Q 65 161 64 176 L 64 176 Z M 80 224 Q 65 225 64 240 Q 65 255 80 256 L 144 256 L 144 256 Q 159 255 160 240 Q 159 225 144 224 L 80 224 L 80 224 Z" />
                    </svg>
                  </Button>
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
                            ) : doc.stage === "Postprocessing" ? (
                              <div className="app-loader" />
                            ) : (
                              <>
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
                              </>
                            )}
                            <Icons.Document />
                          </Tab>
                        ))}
                      </TabList>
                      <Button
                        className={"react-aria-Button new-document"}
                        onPress={newDocument}
                        style={{
                          left: newTabPositionRect?.right || 0,
                          top: newTabPositionRect?.top,
                          opacity: hideNewTab ? "0" : "1",
                        }}
                      >
                        <svg
                          height="16"
                          xmlns="http://www.w3.org/2000/svg"
                          aria-hidden="true"
                          viewBox="0 0 448 512"
                        >
                          <path
                            d="M256 80c0-17.7-14.3-32-32-32s-32 14.3-32 32v144H48c-17.7 0-32 14.3-32 32s14.3 32 32 32h144v144c0 17.7 14.3 32 32 32s32-14.3 32-32V288h144c17.7 0 32-14.3 32-32s-14.3-32-32-32H256V80z"
                            fill="currentColor"
                          />
                        </svg>
                      </Button>
                    </Tabs>
                  ) : (
                    <Select
                      wide
                      selectedKey={appState.activeDocumentId}
                      onSelectionChange={(value) => setActiveDocumentId(value)}
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
                    <svg
                      xmlns="http://www.w3.org/2000/svg"
                      aria-hidden="true"
                      viewBox="0 0 24 24"
                      height="12"
                    >
                      <g
                        style={{
                          fill: "none",
                          stroke: "currentcolor",
                          strokeWidth: 2,
                          strokeLinecap: "round",
                          strokeLinejoin: "round",
                          strokeMiterlimit: 10,
                        }}
                      >
                        <path d="M15 3h6v6M21 3 11 13M19 18v1a2 2 0 0 1-2 2H5a2 2 0 0 1-2-2V7a2 2 0 0 1 2-2h1" />
                      </g>
                    </svg>
                    <span>Export</span>
                  </Button>
                </div>
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
              {activeDocument ? <ExportModal /> : null}
            </div>
          ) : (
            <div className="app-overlay app-loader" />
          )}
        </HeatmapRefProvider>
      </AppStateContext.Provider>
    </ErrorBoundary>
  );
};
