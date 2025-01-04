import React from "react";
import {
  type AppState,
  type DocState,
  type SetAppState,
  docStateSchema,
  findDoc,
  initialDocState,
} from "../appState";
import { partialSafeParse } from "../zodUtils";

export const useDocState = (
  docId: string,
  appState: AppState,
  setAppState: SetAppState,
) => {
  const docState = React.useMemo(() => {
    const state = findDoc(docId, appState.documents);

    if (!state) {
      throw new Error(`Expected document to be found: ${docId}`);
    }

    if (state.parsed) {
      return state;
    }

    try {
      const parsedState = partialSafeParse(docStateSchema, state);
      const validData = parsedState.validData;
      const merged: DocState = {
        ...initialDocState,
        ...validData,
        distribution: {
          ...initialDocState.distribution,
          ...validData.distribution,
          histogram: {
            ...initialDocState.distribution.histogram,
            ...validData.distribution?.histogram,
          },
          violin: {
            ...initialDocState.distribution.violin,
            ...validData.distribution?.violin,
          },
          raincloud: {
            ...initialDocState.distribution.raincloud,
            ...validData.distribution?.raincloud,
          },
        },
        heatmap: {
          ...initialDocState.heatmap,
          ...validData.heatmap,
        },
      };
      if ("compute_stats" in state) {
        merged.compute_stats = state.compute_stats;
      }
      if (parsedState.error) {
        console.warn(parsedState.error);
      }
      setAppState((prev) => ({
        ...prev,
        documents: prev.documents.map((doc) =>
          doc.id === docId ? { ...prev, ...merged, parsed: true } : doc,
        ),
      }));
      return merged;
    } catch (e) {
      return initialDocState;
    }
  }, [docId, appState.documents, setAppState]);

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
