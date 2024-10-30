// HistogramPlot.tsx
import React from "react";
import { PlotData, Layout } from "plotly.js-dist-min";
import Plotly from "plotly.js-dist-min";
import createPlotlyComponent from "react-plotly.js/factory";

const Plot = createPlotlyComponent(Plotly);

interface HistogramPlotProps {
  x: number[];
  y: number[];
  settings: {
    color: string;
    outlineColor: string;
    outlineWidth: number;
  };
}

const HistogramPlot: React.FC<HistogramPlotProps> = ({ x, y, settings }) => {
  const trace: Partial<PlotData> = {
    type: "bar",
    x: x,
    y: y,
    marker: {
      color: settings.color,
      line: {
        width: settings.outlineWidth,
        color: settings.outlineColor,
      },
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
    displayModeBar: false, // Hides the mode bar
    scrollZoom: true,
    displaylogo: false, // Hides the Plotly logo
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
