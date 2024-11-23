import React from "react";
import { Histogram } from "./Histogram";
import { Violin } from "./Violin";
import { Raincloud } from "./Raincloud";
import { DistributionData } from "../plotTypes";

export const Distribution = ({
  data,
  tickText,
}: {
  data: DistributionData | undefined;
  tickText: string[];
  footer?: React.ReactNode;
}) => {
  if (!data) {
    return (
      <div className="app-main centered">
        <h1>No data available</h1>
      </div>
    );
  }
  const [view, setView] = React.useState<"histogram" | "violin" | "raincloud">(
    "violin",
  ); // Switch the default here until we have a switcher

  return (
    <>
      {view === "histogram" && <Histogram data={data} />}
      {view === "violin" && <Violin data={data} tickText={tickText} />}
      {view === "raincloud" && <Raincloud data={data} />}
    </>
  );
};
