import React from "react";
import Plotly from "plotly.js-dist-min";
import { Layout, PlotData } from "plotly.js-dist-min";
import createPlotlyComponent from "react-plotly.js/factory";
import {DistributionData,DistributionSettings } from "../plotTypes";


const Plot = createPlotlyComponent(Plotly);

export interface HistogramPlotSettings {
  color: string;
  outlineColor: string;
  outlineWidth: number;
}



const HistogramPlot: React.FC<DistributionData> = ({ x, y, gc_stats, len_stats, raw_mat, flat_mat }) => {
  const trace: Partial<PlotData> = {
    type: "histogram",
    x: raw_mat,
    marker: {
      color: "pink",
    },
    name: "Histogram",
    hovertemplate: "Percent Identity: %{x}<br>Proportion: %{y}<extra></extra>",
  };

  const layout: Partial<Layout> = {
    margin: {
      l: 0,
      r: 0,
      t: 0,
      b: 0,
    },
    paper_bgcolor: "transparent",
    plot_bgcolor: "transparent",
    xaxis: {
      showline: false,
      zeroline: false,
      showgrid: false,
      showticklabels: false,
    },
    yaxis: {
      showline: false,
      zeroline: false,
      showgrid: false,
      showticklabels: false,
    },
  };

  const config = {
    responsive: true,
    displayModeBar: false, 
    scrollZoom: true,
    displaylogo: false,
  };

  return (
    <div className="plot-container">
      <Plot
        data={[trace]}
        layout={layout}
        config={config}
        style={{ width: "100%", height: "100%" }}
      />
    </div>
  );
};

export default HistogramPlot;
