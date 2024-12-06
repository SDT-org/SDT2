import React from "react";
import {
  type AppState,
  AppStateContext,
  clientStateKey,
  initialAppState,
} from "../appState";
import { useAppBlur } from "../hooks/appBlur";
import { useWaitForPywebview } from "../hooks/usePywebviewReadyEvent";
import { useSaveState } from "../hooks/useSaveState";
import { useShortcutKeys } from "../hooks/useShortcutKeys";
import { useStartProcessData } from "../hooks/useStartProcessData";
import { useSyncState } from "../hooks/useSyncState";
import { ErrorBoundary } from "./ErrorBoundary";
import { ExportModal } from "./ExportModal";
import { Loader } from "./Loader";
import { MainMenu } from "./Menu";
import { Runner } from "./Runner";
import { Viewer } from "./Viewer";

export const App = () => {
  const restoreInitialState = React.useCallback(() => {
    setLoading(true);
    const savedClient = localStorage.getItem(clientStateKey);

    window.pywebview.api.get_state().then((data) => {
      setAppState((prev) => ({
        ...data,
        client: {
          ...prev.client,
          ...(savedClient ? JSON.parse(savedClient) : {}),
          showExportModal: false,
        },
      }));
      setLoading(false);
      setInitialized(true);
    });
  }, []);

  const [appState, setAppState] = React.useState<AppState>(initialAppState);
  const [loading, setLoading] = React.useState(true);
  const [initialized, setInitialized] = React.useState(false);
  const startProcessData = useStartProcessData(appState);
  useSyncState(setAppState);
  useShortcutKeys(appState, setAppState);
  useAppBlur();
  useWaitForPywebview(restoreInitialState);
  useSaveState(initialized, appState);

  const commonViewProps = {
    appState,
    setAppState,
    mainMenu: <MainMenu />,
  };

  const APP_VIEWS: { [K in AppState["view"]]: React.ReactElement } = {
    runner: (
      <Runner
        {...commonViewProps}
        startProcessData={startProcessData}
        appState={appState}
      />
    ),
    loader: <Loader {...commonViewProps} />,
    viewer: <Viewer {...commonViewProps} />,
  };

  return (
    <ErrorBoundary appState={appState} setAppState={setAppState}>
      <AppStateContext.Provider value={{ appState, setAppState }}>
        {loading ? <div className="app-overlay app-loader" /> : null}
        {initialized && !loading ? APP_VIEWS[appState?.view || "viewer"] : null}
        <ExportModal />
      </AppStateContext.Provider>
    </ErrorBoundary>
  );
};
