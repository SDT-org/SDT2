import * as React from "react";
import {
  AppState,
  AppStateContext,
  initialAppState,
  syncAppState,
} from "../appState";
import { Runner } from "./Runner";
import { Loader } from "./Loader";
import { Viewer } from "./Viewer";
import { ErrorBoundary } from "./ErrorBoundary";
import { MainMenu, MainMenuProps } from "./Menu";
import { ExportModal } from "./ExportModal";

export const App = () => {
  const [appState, setAppState] = React.useState<AppState>(initialAppState);
  const mainMenuCallbacks: MainMenuProps = {
    appState,
    onNew: () => {
      if (
        appState.view !== "runner" &&
        confirm("Are you sure? Current results will be cleared.")
      ) {
        window.pywebview.api.reset_state().then(() => {
          // Reset client state as well
          setAppState((previous) => {
            return {
              ...initialAppState,
              client: { ...previous.client },
            };
          });
        });
      }
    },
    onExport: () =>
      setAppState((previous) => {
        return {
          ...previous,
          client: { ...previous.client, showExportModal: true },
        };
      }),
    onSettings: () => {},
    onAbout: window.pywebview.api.show_about,
    onManual: window.pywebview.api.show_manual,
    onExit: () => {
      if (confirm("Are you sure you want to exit?")) {
        window.pywebview.api.close_app();
      }
    },
  };
  const commonViewProps = {
    appState,
    setAppState,
    mainMenu: <MainMenu {...mainMenuCallbacks} />,
  };
  const [showDebugState, setShowDebugState] = React.useState(false);

  const APP_VIEWS: { [K in AppState["view"]]: React.ReactElement } = {
    runner: <Runner {...commonViewProps} />,
    loader: <Loader {...commonViewProps} />,
    viewer: <Viewer {...commonViewProps} />,
  };

  if (typeof window.syncAppState !== "function") {
    window.syncAppState = () => syncAppState(setAppState);
    window.syncAppState();
  }

  if (appState.debug) {
    document.addEventListener("keydown", (event) => {
      if (event.key === "d") {
        setShowDebugState(!showDebugState);
      }
    });
  }

  React.useEffect(() => {
    const handleKeydown = (event: KeyboardEvent) => {
      if ((event.ctrlKey || event.metaKey) && event.key === "o") {
        event.preventDefault();
        window.pywebview.api.open_file_dialog();
      }
    };
    document.addEventListener("keydown", handleKeydown);

    return () => {
      document.removeEventListener("keydown", handleKeydown);
    };
  }, []);

  return (
    <ErrorBoundary appState={appState} setAppState={setAppState}>
      <AppStateContext.Provider value={{ appState, setAppState }}>
        {APP_VIEWS[appState?.view || "viewer"]}
        <ExportModal />
        {showDebugState ? (
          <pre>AppState{JSON.stringify(appState, null, 2)}</pre>
        ) : null}
      </AppStateContext.Provider>
    </ErrorBoundary>
  );
};
