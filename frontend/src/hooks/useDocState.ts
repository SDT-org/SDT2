import React from "react";
import {
  type AppState,
  type DocState,
  type SetAppState,
  findDoc,
  initialDocState,
} from "../appState";

export const useDocState = (
  docId: string,
  appState: AppState,
  setAppState: SetAppState,
) => {
  const foundState = findDoc(docId, appState);

  if (!foundState) {
    throw new Error(
      `Expected document to be found: ${appState.activeDocumentId}`,
    );
  }

  const docState = { ...initialDocState, ...foundState };

  const setDocState = React.useCallback(
    (nextDoc: (prevDoc: DocState) => DocState) =>
      setAppState((prev) => ({
        ...prev,
        documents: prev.documents.map((prevDoc) =>
          prevDoc.filename[0] === docId ? nextDoc(prevDoc) : prevDoc,
        ),
      })),
    [docId, setAppState],
  );

  const updateDocState = React.useCallback(
    (values: Partial<DocState>) => setDocState((prev) => ({ ...prev, values })),
    [setDocState],
  );

  return { docState, setDocState, updateDocState };
};
