import React from "react";
import { Histogram } from "./Histogram";
import { Violin } from "./Violin";
import { Raincloud } from "./Raincloud";
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
    "histogram" | "violin" | "raincloud" 
  >("raincloud"); // Switch the default here until we have a switcher

  return (
    <>
      {view === "histogram" && <Histogram data={data} />}
      {view === "violin" && <Violin data={data} />}
      {view === "raincloud" && <Raincloud data={data} />}
    </>
  );
};
