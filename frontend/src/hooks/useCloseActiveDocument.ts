import type { AppState, SetAppState } from "../appState";
import { useCloseDocument } from "./useCloseDocument";

export const useCloseActiveDocument = (
  appState: AppState,
  setAppState: SetAppState,
) => {
  const closeDocument = useCloseDocument(appState, setAppState);

  return () => {
    const multipleDocuments = appState.documents.length > 1;

    closeDocument(appState.activeDocumentId);
    return multipleDocuments;
  };
};
