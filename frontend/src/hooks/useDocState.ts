import React from "react";
import {
  type AppState,
  type DocState,
  type SetAppState,
  findDoc,
} from "../appState";

export const useDocState = (
  docId: string,
  appState: AppState,
  setAppState: SetAppState,
) => {
  const docState = React.useMemo(() => {
    const foundDoc = findDoc(docId, appState.documents);

    if (!foundDoc) {
      throw new Error(`Expected document to be found: ${docId}`);
    }

    return foundDoc;
  }, [docId, appState.documents]);

  const setDocState = React.useCallback(
    (nextDoc: (prevDoc: DocState) => DocState, markModified = true) => {
      const isValidType = [docState.filetype].includes("application/vnd.sdt");
      const extraValues = isValidType && markModified ? { modified: true } : {};

      return setAppState((prev) => ({
        ...prev,
        documents: prev.documents.map((prevDoc) =>
          prevDoc.id === docId
            ? {
                ...nextDoc(prevDoc),
                ...extraValues,
              }
            : prevDoc,
        ),
      }));
    },
    [docId, setAppState, docState.filetype],
  );

  const updateDocState = React.useCallback(
    (values: Partial<DocState>, markModified = true) => {
      setDocState((prev) => ({ ...prev, ...values }), markModified);
    },
    [setDocState],
  );

  return { docState, setDocState, updateDocState };
};
