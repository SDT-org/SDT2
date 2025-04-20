// import React from "react";
// import type { ExportImageData } from "../appState";

export const useExportData = () => {
  // const [data, setData] = React.useState<ExportImageData>();
  // return [data, setData] as const;
};

// maybe this should be a thing that is in a context that provides "is exporting" flag and a function to call when a component renders (in an effect, or in heatmap, when onRenderComplete) to just send the image data to the backend each time. this is likely less memory intensive than making multiple strings of data to send at once?
