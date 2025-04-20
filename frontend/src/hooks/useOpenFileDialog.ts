import React from "react";
import type { AppState, SetAppState } from "../appState";
import { services } from "../services";

const useOpenFileDialog = (appState: AppState, setAppState: SetAppState) => {
  return React.useCallback(
    (docId?: string) => {
      services.openFileDialog(docId, appState.lastDataFilePath).then((data) => {
        const [success, result] = data;
        if (!success) {
          return;
        }
        setAppState((prev) => ({
          ...prev,
          activeDocumentId: result[0],
          lastDataFilePath: result[1],
        }));
      });
    },
    [appState.lastDataFilePath, setAppState],
  );
};

export default useOpenFileDialog;
