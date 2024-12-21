import React from "react";
import { type AppState, AppStateContext, initialAppState } from "../appState";
import { useAppBlur } from "../hooks/appBlur";
import { useWaitForPywebview } from "../hooks/usePywebviewReadyEvent";
import { useSaveState } from "../hooks/useSaveState";
import { useShortcutKeys } from "../hooks/useShortcutKeys";
import { useStartRun } from "../hooks/useStartProcessData";
import { useSyncState } from "../hooks/useSyncState";
import { restoreClientState } from "../restoreClientState";
import { ErrorBoundary } from "./ErrorBoundary";
import { ExportModal } from "./ExportModal";
import { Loader } from "./Loader";
import { MainMenu } from "./Menu";
import { Runner } from "./Runner";
import { Viewer } from "./Viewer";

export const App = () => {
  const restoreInitialState = React.useCallback(() => {
    setLoading(true);

    window.pywebview.api.get_state().then((data) => {
      setAppState((prev) => ({
        ...data,
        client: restoreClientState(prev.client),
        showExportModal: false,
      }));
      setLoading(false);
      setInitialized(true);
    });
  }, []);

  const [appState, setAppState] = React.useState<AppState>(initialAppState);
  const [loading, setLoading] = React.useState(true);
  const [initialized, setInitialized] = React.useState(false);
  const startRun = useStartRun(appState);
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
      <Runner {...commonViewProps} startRun={startRun} appState={appState} />
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
