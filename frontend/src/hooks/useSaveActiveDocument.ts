import { type AppState, type SetAppState, findDoc } from "../appState";
import { isSDTFile } from "../helpers";
import { saveDocument } from "../services/documents";

export const useSaveActiveDocument = (
  appState: AppState,
  setAppState: SetAppState,
) => {
  return (saveAs = false) => {
    const doc = findDoc(appState.activeDocumentId, appState.documents);
    const saveable = doc && doc.view === "viewer";

    if (saveable) {
      saveDocument(doc, saveAs || !isSDTFile(doc.filetype)).then(() => {
        setAppState((prev) => ({
          ...prev,
          documents: prev.documents.map((prevDoc) => ({
            ...prevDoc,
            ...(prevDoc.id === doc.id ? { modified: false } : {}),
          })),
        }));
      });
    }

    return saveable;
  };
};
