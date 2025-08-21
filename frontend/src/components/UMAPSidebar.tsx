import type React from "react";
import type { UMAPSettings } from "../plotTypes";
import { Slider } from "./primitives/Slider";

interface UMAPSidebarProps {
  settings: UMAPSettings;
  updateSettings: (values: Partial<UMAPSettings>) => void;
  leftSidebarCollapsed: boolean;
}

export const UMAPSidebar: React.FC<UMAPSidebarProps> = ({
  settings,
  updateSettings,
}) => {
  const handleParameterChange = (param: string, value: number | boolean) => {
    updateSettings({ [param]: value });
  };

  return (
    <div className="app-sidebar app-sidebar-right umap-sidebar">
      <div className="app-sidebar-body">
        <div className="sidebar-section">
          <h3>UMAP Parameters</h3>
          <p className="sidebar-description">
            UMAP (Uniform Manifold Approximation and Projection) reduces
            high-dimensional distance data to 2D for visualization.
          </p>

          <div className="sidebar-item">
            <label htmlFor="n_neighbors">n_neighbors</label>
            <p className="param-description">
              Controls how UMAP balances local vs global structure. Lower values
              focus on local structure.
            </p>
            <Slider
              id="n_neighbors"
              value={settings.n_neighbors}
              onChangeEnd={(value: number) =>
                handleParameterChange("n_neighbors", value)
              }
              minValue={2}
              maxValue={200}
              step={1}
            />
            <span className="value-display">{settings.n_neighbors}</span>
          </div>

          <div className="sidebar-item">
            <label htmlFor="min_dist">min_dist</label>
            <p className="param-description">
              Minimum distance between points. Lower values create tighter
              clusters.
            </p>
            <Slider
              id="min_dist"
              value={settings.min_dist}
              onChangeEnd={(value: number) =>
                handleParameterChange("min_dist", value)
              }
              minValue={0}
              maxValue={1}
              step={0.01}
            />
            <span className="value-display">
              {settings.min_dist.toFixed(2)}
            </span>
          </div>
        </div>

        <div className="sidebar-section">
          <h3>HDBSCAN Clustering</h3>
          <p className="sidebar-description">
            HDBSCAN finds clusters of varying densities in the data. Points not
            belonging to any cluster are marked as noise.
          </p>

          <div className="sidebar-item">
            <label htmlFor="minClusterSize">Min Cluster Size</label>
            <p className="param-description">
              Minimum number of points required to form a cluster. Smaller
              values find more clusters.
            </p>
            <Slider
              id="minClusterSize"
              value={settings.minClusterSize}
              onChangeEnd={(value: number) =>
                handleParameterChange("minClusterSize", value)
              }
              minValue={2}
              maxValue={50}
              step={1}
            />
            <span className="value-display">{settings.minClusterSize}</span>
          </div>

          <div className="sidebar-item">
            <label htmlFor="clusterEpsilon">
              Cluster Selection Epsilon (%)
            </label>
            <p className="param-description">
              Minimum similarity percentage for clustering. Points must be at
              least this similar to belong to the same cluster. 0 uses the full
              hierarchy.
            </p>
            <Slider
              id="clusterEpsilon"
              value={settings.clusterEpsilon}
              onChangeEnd={(value: number) =>
                handleParameterChange("clusterEpsilon", value)
              }
              minValue={0}
              maxValue={100}
              step={0.5}
            />
            <span className="value-display">
              {settings.clusterEpsilon.toFixed(1)}
            </span>
          </div>
        </div>
      </div>
    </div>
  );
};
