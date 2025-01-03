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
    (nextDoc: (prevDoc: DocState) => DocState) =>
      setAppState((prev) => ({
        ...prev,
        documents: prev.documents.map((prevDoc) =>
          prevDoc.id === docId
            ? { ...nextDoc(prevDoc), modified: true }
            : prevDoc,
        ),
      })),
    [docId, setAppState],
  );

  const updateDocState = React.useCallback(
    (values: Partial<DocState>) =>
      setDocState((prev) => ({ ...prev, ...values })),
    [setDocState],
  );

  return { docState, setDocState, updateDocState };
};
