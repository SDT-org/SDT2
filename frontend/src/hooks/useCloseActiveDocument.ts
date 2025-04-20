import type { AppState, SetAppState } from "../appState";
import { useCloseDocument } from "./useCloseDocument";

export const useCloseActiveDocument = (
  appState: AppState,
  setAppState: SetAppState,
) => {
  const closeDocument = useCloseDocument(appState, setAppState);

  return async () => {
    return await closeDocument(appState.activeDocumentId);
  };
};
