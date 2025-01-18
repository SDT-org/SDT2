import React from "react";
import { HeatmapRefContext, type HeatmapRefType } from "../appState";

export const HeatmapRefProvider = ({
  children,
}: { children: React.ReactNode }) => {
  const heatmapRef = React.useRef<HeatmapRefType | null>(null);
  return (
    <HeatmapRefContext.Provider value={heatmapRef}>
      {children}
    </HeatmapRefContext.Provider>
  );
};
