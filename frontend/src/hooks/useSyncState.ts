import React from "react";
import {
  type SetAppState,
  type SyncProgressEvent,
  type SyncStateEvent,
  findDoc,
} from "../appState";

export const useSyncState = (setAppState: SetAppState) => {
  React.useEffect(() => {
    const handler = (event: SyncStateEvent) => {
      console.count("sync-state");
      setAppState((prev) => ({
        ...prev,
        ...event.detail.state,
        documents: event.detail.state.documents.map((beDoc) => {
          const feDoc = findDoc(beDoc.id, prev.documents);
          return {
            ...beDoc,
            dataView: feDoc?.dataView || beDoc.dataView,
            heatmap: {
              ...beDoc.heatmap,
              ...feDoc?.heatmap,
            },
            distribution: {
              ...beDoc.distribution,
              ...feDoc?.distribution,
            },
          };
        }),
      }));
    };
    document.addEventListener("sync-state", handler);
    return () => document.removeEventListener("sync-state", handler);
  }, [setAppState]);

  React.useEffect(() => {
    const handler = (event: SyncProgressEvent) => {
      setAppState((prev) => ({
        ...prev,
        documents: prev.documents.map((doc) => ({
          ...doc,
          ...(event.detail.id === doc.id ? event.detail : null),
        })),
      }));
    };

    document.addEventListener("sync-progress", handler);
    return () => document.removeEventListener("sync-progress", handler);
  }, [setAppState]);
};
