import React from "react";
import { Histogram } from "./Histogram";
import { Violin } from "./Violin";
import { DistributionData } from "../plotTypes";

export const Distribution = ({
  data,
}: {
  data: DistributionData | undefined;
  footer?: React.ReactNode;
}) => {
  if (!data) {
    return (
      <div className="app-main centered">
        <h1>No data available</h1>
      </div>
    );
  }
  const [view, setView] = React.useState<
    "histogram" | "line" | "violin" | "scatter" | "box"
  >("histogram"); // Switch the default here until we have a switcher

  return (
    <>
      {view === "histogram" && <Histogram data={data} />}
      {view === "line" && <Histogram data={data} />}
      {view === "violin" && <Violin data={data} />}
      {view === "scatter" && <Histogram data={data} />}
    </>
  );
};
