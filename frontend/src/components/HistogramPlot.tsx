import React from "react";
import Plotly from "plotly.js-dist-min";
import { Layout, PlotData } from "plotly.js-dist-min";
import createPlotlyComponent from "react-plotly.js/factory";
import {DistributionData} from "../plotTypes";

// need to import settings

const Plot = createPlotlyComponent(Plotly);

export interface HistogramPlotSettings {
  color: string;
  outlineColor: string;
  outlineWidth: number;
}



const HistogramPlot = ({ x, y, gc_stats, len_stats, raw_mat, }: DistributionData) => {
  const trace: Partial<PlotData> = {
    x: raw_mat,
    type: "histogram",
        histnorm: 'probability',
        marker: {
          color: 'pink',
        },
    name: "Histogram",
    hovertemplate: "Percent Identity: %{x}<br>Proportion: %{y}<extra></extra>",
  };

  console.log('raw_mat', raw_mat);
  console.log('len_stats', len_stats);
  console.log('gc_stats', gc_stats);

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
