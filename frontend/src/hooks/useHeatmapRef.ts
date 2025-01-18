import React from "react";
import { HeatmapRefContext } from "../appState";

export const useHeatmapRef = () => {
  const context = React.useContext(HeatmapRefContext);
  if (!context) {
    throw new Error("useHeatmapRef must be used within a HeatmapRefProvider");
  }
  return context;
};
