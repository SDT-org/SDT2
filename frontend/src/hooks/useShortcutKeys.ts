import React from "react";
import { type AppState, type SetAppState, findDoc } from "../appState";
import { isSDTFile } from "../helpers";
import { saveDocument } from "../services/documents";
import useOpenFileDialog from "./useOpenFileDialog";

export const useShortcutKeys = (
  appState: AppState,
  setAppState: SetAppState,
) => {
  const openFileDialog = useOpenFileDialog(appState, setAppState);

  React.useEffect(() => {
    const handleKeydown = (event: KeyboardEvent) => {
      if (event.ctrlKey || event.metaKey) {
        if (event.key === "o") {
          event.preventDefault();
          openFileDialog();
        }
        if (event.key === "s") {
          event.preventDefault();
          const doc = findDoc(appState.activeDocumentId, appState.documents);
          if (doc && doc.view === "viewer") {
            saveDocument(doc, !isSDTFile(doc.filetype));
          }
        }
      }
    };

    document.addEventListener("keydown", handleKeydown);

    return () => document.removeEventListener("keydown", handleKeydown);
  }, [openFileDialog, appState.activeDocumentId, appState.documents]);
};
