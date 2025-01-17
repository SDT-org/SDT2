import React from "react";
import type { AppState, SetAppState } from "../appState";
import { useCloseActiveDocument } from "./useCloseActiveDocument";
import useNewDocument from "./useNewDocument";
import useOpenFileDialog from "./useOpenFileDialog";
import { useSaveActiveDocument } from "./useSaveActiveDocument";

export const useShortcutKeys = (
  appState: AppState,
  setAppState: SetAppState,
) => {
  const openFileDialog = useOpenFileDialog(appState, setAppState);
  const newDocument = useNewDocument(setAppState);
  const closeActiveDocument = useCloseActiveDocument(appState, setAppState);
  const saveActiveDocument = useSaveActiveDocument(appState, setAppState);

  const keyMap = React.useMemo(
    () => ({
      o: openFileDialog,
      n: newDocument,
      s: saveActiveDocument,
      w: closeActiveDocument,
    }),
    [openFileDialog, newDocument, closeActiveDocument, saveActiveDocument],
  );

  React.useEffect(() => {
    const handleKeydown = (event: KeyboardEvent) => {
      if (
        (event.ctrlKey || event.metaKey) &&
        Object.keys(keyMap).includes(event.key)
      ) {
        const fn = keyMap[event.key as keyof typeof keyMap];
        if (fn?.()) {
          event.preventDefault();
        }
      }
    };

    document.addEventListener("keydown", handleKeydown);

    return () => document.removeEventListener("keydown", handleKeydown);
  }, [keyMap]);
};
