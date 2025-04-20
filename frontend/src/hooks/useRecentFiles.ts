import React from "react";
import type { AppState, DocState, SetAppState } from "../appState";
import { openFile } from "../services/files";

export const useRecentFiles = (
  appState: AppState,
  setAppState: SetAppState,
) => {
  const fetchRecentFiles = React.useCallback(() => {
    window.pywebview.api.app_settings().then((data) => {
      setAppState((prev) => ({ ...prev, recentFiles: data.recent_files }));
    });
  }, [setAppState]);

  React.useEffect(() => {
    fetchRecentFiles();
  }, [fetchRecentFiles]);

  const openRecentFile = React.useCallback(
    async (file: string, doc: DocState) => {
      const isNewDocument = !doc.filetype && !doc.basename;
      const existingDocument = appState.documents.find(
        (document) => document.filename === file,
      );

      if (existingDocument) {
        setAppState((prev) => ({
          ...prev,
          activeDocumentId: existingDocument.id,
        }));
        return;
      }

      const docId = isNewDocument ? doc.id : undefined;
      await openFile(file, docId);
      fetchRecentFiles();
    },
    [appState.documents, setAppState, fetchRecentFiles],
  );

  return openRecentFile;
};
