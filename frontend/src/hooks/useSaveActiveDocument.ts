import { type AppState, findDoc } from "../appState";
import { isSDTFile } from "../helpers";
import { saveDocument } from "../services/documents";

export const useSaveActiveDocument = (appState: AppState) => {
  return () => {
    const doc = findDoc(appState.activeDocumentId, appState.documents);
    const saveable = doc && doc.view === "viewer";

    if (saveable) {
      saveDocument(doc, !isSDTFile(doc.filetype));
    }

    return saveable;
  };
};
