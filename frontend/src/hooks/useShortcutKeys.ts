import React from "react";
import type { AppState, SetAppState } from "../appState";
import useOpenFileDialog from "./useOpenFileDialog";

export const useShortcutKeys = (
  appState: AppState,
  setAppState: SetAppState,
) => {
  const openFileDialog = useOpenFileDialog(appState, setAppState);

  React.useEffect(() => {
    const handleKeydown = (event: KeyboardEvent) => {
      if ((event.ctrlKey || event.metaKey) && event.key === "o") {
        event.preventDefault();
        openFileDialog();
      }
    };
    document.addEventListener("keydown", handleKeydown);

    return () => document.removeEventListener("keydown", handleKeydown);
  }, [openFileDialog]);
};
