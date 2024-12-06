import React from "react";
import type { SetAppState, SyncStateEvent } from "../appState";

export const useSyncState = (setAppState: SetAppState) =>
  React.useEffect(() => {
    const handler = (event: SyncStateEvent) => {
      setAppState((prev) => ({ ...prev, ...event.detail.state }));
    };
    document.addEventListener("sync-state", handler);
    return () => document.removeEventListener("sync-state", handler);
  }, [setAppState]);
