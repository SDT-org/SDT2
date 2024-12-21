import React from "react";
import type { AppState } from "../appState";
import { services } from "../services";

export const useStartRun = (appState: AppState) =>
  React.useCallback(() => {
    services.startRun(appState);
  }, [appState]);
