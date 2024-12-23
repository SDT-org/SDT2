import React from "react";
import type { AppState, DocState } from "../appState";
import { services } from "../services";

export const useStartRun = (docState: DocState, appState: AppState) =>
  React.useCallback(() => {
    services.startRun(docState.id, appState);
  }, [docState, appState]);
