import React from "react";
import type { AppState, SetAppState } from "../appState";

export const useShortcutKeys = (appState: AppState, setAppState: SetAppState) =>
  React.useEffect(() => {
    const handleKeydown = (event: KeyboardEvent) => {
      if ((event.ctrlKey || event.metaKey) && event.key === "o") {
        event.preventDefault();
        window.pywebview.api
          .open_file_dialog(appState.client.lastDataFilePath)
          .then((data) =>
            setAppState((prev) => ({
              ...prev,
              client: { ...prev.client, lastDataFilePath: data },
            })),
          );
      }
    };
    document.addEventListener("keydown", handleKeydown);

    return () => document.removeEventListener("keydown", handleKeydown);
  }, [appState, setAppState]);
