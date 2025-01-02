import type { AppState, DocState, SetAppState } from "../appState";
import { services } from "../services";
import useNewDocument from "./useNewDocument";

export const useCloseDocument = (
  appState: AppState,
  setAppState: SetAppState,
) => {
  return (docId: string) => {
    const newDocument = useNewDocument(setAppState);

    if ([0, 1].includes(appState.documents.length)) {
      return newDocument().then(() => services.closeDocument(docId));
    }

    const findDocIndex = (docs: DocState[]) =>
      docs.findIndex((doc) => doc.id === docId);

    if (appState.activeDocumentId === docId) {
      setAppState((prev) => ({
        ...prev,
        activeDocumentId:
          prev.documents[findDocIndex(prev.documents) + 1]?.id ||
          prev.documents[findDocIndex(prev.documents) - 1]?.id ||
          "",
      }));
    }

    return services.closeDocument(docId);
  };
};
