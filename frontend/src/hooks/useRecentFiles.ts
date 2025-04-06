import React from "react";
import type { DocState, SetAppState } from "../appState";
import { openFile } from "../services/files";

export const useRecentFiles = (setAppState: SetAppState) => {
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
      await openFile(file, doc.id);
      fetchRecentFiles();
    },
    [fetchRecentFiles],
  );

  return openRecentFile;
};
