import * as React from "react";
import { AppState, AppStateContext, initialAppState } from "../appState";
import { Runner } from "./Runner";
import { Loader } from "./Loader";
import { Viewer } from "./Viewer";
import { ErrorBoundary } from "./ErrorBoundary";
import { MainMenu, MainMenuProps } from "./Menu";
import { ExportModal } from "./ExportModal";

export const App = () => {
  const [appState, setAppState] = React.useState<AppState>(initialAppState);
  const [debugState, setDebugState] = React.useState("");
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
    onOpen: () => {
      if (
        appState.view !== "runner" &&
        !confirm("Are you sure? Current results will be cleared.")
      ) {
        return;
      }

      window.pywebview.api.open_file_dialog();
    },
    onExport: () =>
      setAppState((previous) => {
        return {
          ...previous,
          client: { ...previous.client, showExportModal: true },
        };
      }),
    onAbout: () => window.pywebview.api.show_about(),
    onManual: () => window.pywebview.api.show_manual(),
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

  const startProcessData = React.useCallback(() => {
    window.pywebview.api
      .run_process_data({
        cluster_method: appState.client.cluster_method,
        compute_cores: appState.client.compute_cores,
      })
      .catch((e) => {
        if (e.toString().includes("PARASAIL_TRACEBACK")) {
          alert(
            "An error occured while aligning. " +
              "Please ensure you have adequate swap/page size and system memory.\n\n" +
              "Error ID: PARASAIL_TRACEBACK",
          );
          window.pywebview.api.cancel_run();
        } else {
          throw e;
        }
      });
  }, []);

  const APP_VIEWS: { [K in AppState["view"]]: React.ReactElement } = {
    runner: <Runner {...commonViewProps} startProcessData={startProcessData} />,
    loader: <Loader {...commonViewProps} />,
    viewer: <Viewer {...commonViewProps} />,
  };

  if (typeof window.syncAppState !== "function") {
    window.syncAppState = (state: AppState) => {
      setAppState((previous) => {
        return { ...previous, ...state };
      });
    };
  }

  if (appState.debug) {
    document.addEventListener("keydown", (event) => {
      if (event.key === "d") {
        setShowDebugState(!showDebugState);
      }
    });
    (window as any).APP_STATE = appState;
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

  React.useEffect(() => {
    setDebugState(JSON.stringify(appState, null, 2));
  }, [appState]);

  return (
    <ErrorBoundary appState={appState} setAppState={setAppState}>
      <AppStateContext.Provider value={{ appState, setAppState }}>
        {APP_VIEWS[appState?.view || "viewer"]}
        <ExportModal />
        {showDebugState ? <pre>AppState {debugState}</pre> : null}
      </AppStateContext.Provider>
    </ErrorBoundary>
  );
};
