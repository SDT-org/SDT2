import React from "react";
import type { DocState } from "../appState";
import { services } from "../services";

export const useStartRun = (docState: DocState) =>
  React.useCallback(() => services.startRun(docState), [docState]);
