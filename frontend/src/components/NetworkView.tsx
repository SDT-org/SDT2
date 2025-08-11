import Graph from "graphology";
import forceAtlas2 from "graphology-layout-forceatlas2";
import React, { useEffect, useRef, useState } from "react";
import Sigma from "sigma";
import type { DocState, SetDocState } from "../appState";
import { distinctColor } from "../colors";

// Helper function to convert HSL to Hex
function hslToHex(hsl: string): string {
  console.log("Converting HSL to hex:", hsl);
  const match = hsl.match(/hsl\((\d+(?:\.\d+)?),\s*(\d+)%,\s*(\d+)%\)/);
  if (!match || match.length < 4) {
    console.log("Failed to parse HSL, returning default gray");
    return "#cccccc";
  }

  const hue = match[1] || "0";
  const sat = match[2] || "0";
  const light = match[3] || "0";

  const h = Number.parseFloat(hue) / 360;
  const s = Number.parseFloat(sat) / 100;
  const l = Number.parseFloat(light) / 100;

  const hue2rgb = (p: number, q: number, tInput: number) => {
    let t = tInput;
    if (t < 0) t += 1;
    if (t > 1) t -= 1;
    if (t < 1 / 6) return p + (q - p) * 6 * t;
    if (t < 1 / 2) return q;
    if (t < 2 / 3) return p + (q - p) * (2 / 3 - t) * 6;
    return p;
  };

  const q = l < 0.5 ? l * (1 + s) : l + s - l * s;
  const p = 2 * l - q;

  const r = Math.round(hue2rgb(p, q, h + 1 / 3) * 255);
  const g = Math.round(hue2rgb(p, q, h) * 255);
  const b = Math.round(hue2rgb(p, q, h - 1 / 3) * 255);

  const hex = `#${((1 << 24) + (r << 16) + (g << 8) + b).toString(16).slice(1)}`;
  console.log(`Converted ${hsl} to ${hex}`);
  return hex;
}

interface NetworkViewProps {
  docState: DocState;
  setDocState: SetDocState;
  leftSidebarCollapsed: boolean;
}

interface NetworkData {
  nodes: Array<{
    id: string;
    x: number;
    y: number;
    cluster: number;
    degree: number;
  }>;
  edges: Array<{
    source: string;
    target: string;
    weight: number;
  }>;
  stats: {
    total_nodes: number;
    total_edges: number;
    total_clusters: number;
    noise_points: number;
    largest_cluster_size: number;
    smallest_cluster_size: number;
    average_degree: number;
    density: number;
    cluster_sizes: Record<number, number>;
  };
  bounds: {
    x: [number, number];
    y: [number, number];
  };
}

export const NetworkView: React.FC<NetworkViewProps> = ({
  docState,
  setDocState,
  leftSidebarCollapsed: _leftSidebarCollapsed,
}) => {
  const containerRef = useRef<HTMLDivElement>(null);
  const sigmaRef = useRef<Sigma | null>(null);
  const [networkData, setNetworkData] = useState<NetworkData | null>(null);
  const [loading, setLoading] = useState(false);
  const [error, setError] = useState<string | null>(null);
  const [hoveredNode, setHoveredNode] = useState<string | null>(null);
  const [selectedNode, setSelectedNode] = useState<string | null>(null);
  const [isLayoutRunning, setIsLayoutRunning] = useState(false);

  const updateSettings = React.useCallback(
    (values: Partial<DocState["network"]>) =>
      setDocState((prev) => ({
        ...prev,
        network: {
          ...prev.network,
          ...values,
        },
      })),
    [setDocState],
  );

  // Fetch network data
  useEffect(() => {
    if (!docState.id) return;

    const fetchNetworkData = async () => {
      setLoading(true);
      setError(null);

      try {
        const response = await window.pywebview.api.data.get_network_data(
          docState.id,
          {
            similarity_threshold: docState.network.similarityThreshold,
            clustering_method: "louvain",
            layout_method: "spring",
            resolution: docState.network.resolution,
            min_similarity_filter: 85.0, // Increase to 85% to dramatically reduce edges
          },
        );

        console.log("Network data received:", response.data);
        console.log("Nodes:", response.data?.nodes?.length);
        console.log("Edges:", response.data?.edges?.length);
        console.log("Stats:", response.data?.stats);

        setNetworkData(response.data);
      } catch (err) {
        console.error("Error fetching network data:", err);
        setError(
          err instanceof Error ? err.message : "Failed to load network data",
        );
      } finally {
        setLoading(false);
      }
    };

    const timeoutId = setTimeout(fetchNetworkData, 500);
    return () => clearTimeout(timeoutId);
  }, [
    docState.id,
    docState.network.similarityThreshold,
    docState.network.resolution,
  ]);

  // Render network with Sigma.js
  useEffect(() => {
    if (!networkData || !containerRef.current || loading || error) return;

    console.log("Starting graph creation, container:", containerRef.current);
    console.log(
      "NetworkData nodes/edges:",
      networkData.nodes?.length,
      networkData.edges?.length,
    );

    // Skip rendering if there are no edges
    if (networkData.edges.length === 0) {
      console.log("No edges to render, skipping Sigma initialization");
      return;
    }

    // Check container dimensions
    const rect = containerRef.current.getBoundingClientRect();
    console.log("Container dimensions:", rect.width, rect.height);
    if (rect.width === 0 || rect.height === 0) {
      console.error("Container has zero dimensions, cannot render");
      return;
    }

    // Clear previous sigma instance
    if (sigmaRef.current) {
      sigmaRef.current.kill();
      sigmaRef.current = null;
    }

    try {
      // Create graph
      const graph = new Graph();
      console.log("Graph created");

      // Add nodes with proper cluster coloring
      let nodesAdded = 0;
      const clusterColors = new Map<number, string>();

      for (const node of networkData.nodes) {
        // Get or create color for this cluster
        if (!clusterColors.has(node.cluster)) {
          const hslColor =
            node.cluster === 0
              ? "hsl(0, 0%, 80%)"
              : distinctColor(node.cluster);
          const hexColor = hslToHex(hslColor);
          clusterColors.set(node.cluster, hexColor);
          console.log(`Cluster ${node.cluster} assigned color: ${hexColor}`);
        }

        const color = clusterColors.get(node.cluster) || "#cccccc";

        graph.addNode(node.id, {
          x: node.x,
          y: node.y,
          size: Math.max(4, Math.min(15, 4 + Math.log(node.degree + 1) * 2)), // Smaller size range
          label: node.id,
          color: color,
          cluster: node.cluster,
          degree: node.degree,
          borderColor: color,
        });
        nodesAdded++;
      }
      console.log(
        `Added ${nodesAdded} nodes with ${clusterColors.size} different cluster colors`,
      );

      // Add edges
      let edgesAdded = 0;
      for (const edge of networkData.edges) {
        try {
          graph.addEdge(edge.source, edge.target, {
            weight: edge.weight,
            size: Math.max(0.5, Math.min(3, edge.weight / 50)),
            color: "#e0e0e0",
          });
          edgesAdded++;
        } catch (edgeError) {
          console.warn(
            `Failed to add edge ${edge.source} -> ${edge.target}:`,
            edgeError,
          );
        }
      }
      console.log(`Added ${edgesAdded} edges`);

      // Create Sigma instance
      console.log("Creating Sigma instance");
      const sigma = new Sigma(graph, containerRef.current, {
        renderEdgeLabels: false,
        labelDensity: 0.07,
        labelGridCellSize: 60,
        labelRenderedSizeThreshold: 12,
        labelFont: "Arial",
        zoomDuration: 200,
        enableEdgeEvents: true,
        nodeReducer: (node: string, data: Record<string, unknown>) => {
          const res = { ...data };

          if (
            hoveredNode &&
            hoveredNode !== node &&
            graph.neighbors(hoveredNode).indexOf(node) === -1
          ) {
            // biome-ignore lint/complexity/useLiteralKeys: Required for TypeScript index signature
            res["label"] = "";
            // biome-ignore lint/complexity/useLiteralKeys: Required for TypeScript index signature
            res["color"] = "#e0e0e0"; // Lighter gray for non-neighbor nodes
          }
          if (selectedNode === node) {
            // biome-ignore lint/complexity/useLiteralKeys: Required for TypeScript index signature
            res["highlighted"] = true;
            // biome-ignore lint/complexity/useLiteralKeys: Required for TypeScript index signature
            res["size"] =
              (graph.getNodeAttribute(node, "size") as number) * 1.5; // Make selected node bigger
          }
          return res;
        },
        edgeReducer: (edge: string, data: Record<string, unknown>) => {
          const res = { ...data };
          if (hoveredNode && !graph.hasExtremity(edge, hoveredNode)) {
            // biome-ignore lint/complexity/useLiteralKeys: Required for TypeScript index signature
            res["hidden"] = true;
          }
          return res;
        },
      });

      sigmaRef.current = sigma;
      console.log("Sigma instance created successfully");

      // Event listeners
      sigma.on("enterNode", ({ node }: { node: string }) => {
        setHoveredNode(node);
        if (containerRef.current) {
          containerRef.current.style.cursor = "pointer";
        }
      });

      sigma.on("leaveNode", () => {
        setHoveredNode(null);
        if (containerRef.current) {
          containerRef.current.style.cursor = "default";
        }
      });

      sigma.on("clickNode", ({ node }: { node: string }) => {
        setSelectedNode(node === selectedNode ? null : node);
      });

      sigma.on("clickStage", () => {
        setSelectedNode(null);
      });

      // Store sigma instance
      sigmaRef.current = sigma;

      // Cleanup function
      return () => {
        if (sigmaRef.current) {
          sigmaRef.current.kill();
          sigmaRef.current = null;
        }
      };
    } catch (renderError) {
      console.error("Failed to render network:", renderError);
      setError(
        `Failed to render network: ${renderError instanceof Error ? renderError.message : String(renderError)}`,
      );
    }

    // Return cleanup function for all code paths
    return () => {
      if (sigmaRef.current) {
        sigmaRef.current.kill();
        sigmaRef.current = null;
      }
    };
  }, [networkData, loading, error, hoveredNode, selectedNode]);

  // Run ForceAtlas2 layout
  const runForceAtlas2 = () => {
    if (!sigmaRef.current || !networkData || isLayoutRunning) return;

    setIsLayoutRunning(true);
    const graph = sigmaRef.current.getGraph();

    // Run layout
    forceAtlas2.assign(graph, {
      iterations: 100,
      settings: {
        gravity: 1,
        scalingRatio: 10,
        barnesHutOptimize: true,
        barnesHutTheta: 0.5,
        strongGravityMode: false,
        outboundAttractionDistribution: false,
      },
    });

    sigmaRef.current.refresh();
    setIsLayoutRunning(false);
  };

  const resetView = () => {
    if (sigmaRef.current) {
      sigmaRef.current.getCamera().setState({ x: 0.5, y: 0.5, ratio: 1 });
    }
  };

  return (
    <>
      <div className="app-main network-container">
        {loading && (
          <>
            <div className="app-loader app-main-loader" aria-hidden="true" />
            <div className="app-loader-status">
              Loading network visualization...
            </div>
          </>
        )}
        {error && (
          <div
            className="error-message"
            style={{ padding: "20px", color: "red" }}
          >
            {error}
          </div>
        )}
        <div
          ref={containerRef}
          style={{
            width: "100%",
            height: "100%",
            position: "relative",
          }}
        />
        {networkData && networkData.edges.length === 0 && (
          <div
            style={{
              position: "absolute",
              top: "50%",
              left: "50%",
              transform: "translate(-50%, -50%)",
              textAlign: "center",
              padding: "20px",
              background: "rgba(255, 255, 255, 0.9)",
              borderRadius: "8px",
              boxShadow: "0 2px 8px rgba(0,0,0,0.1)",
            }}
          >
            <h3>No connections found</h3>
            <p>
              Your similarity threshold ({docState.network.similarityThreshold}
              %) is too high.
            </p>
            <p>Try lowering the threshold to see more connections.</p>
            <p style={{ fontSize: "12px", color: "#666" }}>
              {networkData.nodes.length} sequences, {networkData.edges.length}{" "}
              edges
            </p>
          </div>
        )}
        {networkData && (
          <div className="network-controls">
            <button
              type="button"
              onClick={runForceAtlas2}
              disabled={isLayoutRunning}
            >
              {isLayoutRunning ? "Running layout..." : "Run ForceAtlas2"}
            </button>
            <button type="button" onClick={resetView}>
              Reset View
            </button>
          </div>
        )}
        {networkData && (
          <div className="network-stats">
            <h3>Louvain Clustering</h3>
            <p>
              <strong>Clusters:</strong> {networkData.stats.total_clusters}
            </p>
            <p>
              <strong>Singletons:</strong> {networkData.stats.noise_points}
            </p>
            <p>
              <strong>Largest cluster:</strong>{" "}
              {networkData.stats.largest_cluster_size} nodes
            </p>
            <p>
              <strong>Edges shown:</strong>{" "}
              {networkData.stats.total_edges.toLocaleString()}
            </p>
            <p>
              <strong>Avg connections:</strong>{" "}
              {networkData.stats.average_degree.toFixed(1)}
            </p>
            {networkData.stats.total_edges > 10000 && (
              <p style={{ color: "orange", fontSize: "12px" }}>
                ⚠️ Many edges - consider increasing threshold
              </p>
            )}
            <div className="cluster-legend" style={{ marginTop: "15px" }}>
              <h4>Cluster Colors</h4>
              <div style={{ maxHeight: "200px", overflowY: "auto" }}>
                {Object.entries(networkData.stats.cluster_sizes)
                  .filter(([cluster]) => cluster !== "0")
                  .sort(([, a], [, b]) => b - a)
                  .slice(0, 20)
                  .map(([cluster, size]) => (
                    <div
                      key={cluster}
                      style={{
                        display: "flex",
                        alignItems: "center",
                        marginBottom: "4px",
                      }}
                    >
                      <div
                        style={{
                          width: "16px",
                          height: "16px",
                          backgroundColor: hslToHex(
                            distinctColor(Number(cluster)),
                          ),
                          marginRight: "8px",
                          borderRadius: "50%",
                        }}
                      />
                      <span style={{ fontSize: "12px" }}>
                        Cluster {cluster}: {size} nodes
                      </span>
                    </div>
                  ))}
                {networkData.stats.noise_points > 0 && (
                  <div
                    style={{
                      display: "flex",
                      alignItems: "center",
                      marginTop: "8px",
                    }}
                  >
                    <div
                      style={{
                        width: "16px",
                        height: "16px",
                        backgroundColor: "#cccccc",
                        marginRight: "8px",
                        borderRadius: "50%",
                      }}
                    />
                    <span style={{ fontSize: "12px" }}>
                      Singletons: {networkData.stats.noise_points} nodes
                    </span>
                  </div>
                )}
              </div>
            </div>
          </div>
        )}
      </div>
      <div className="network-sidebar">
        <h3>Network Settings</h3>
        <div className="control-group">
          <label htmlFor="similarity-threshold">
            Similarity Threshold: {docState.network.similarityThreshold}%
          </label>
          <input
            id="similarity-threshold"
            type="range"
            min="70"
            max="100"
            step="1"
            value={docState.network.similarityThreshold}
            onChange={(e) =>
              updateSettings({ similarityThreshold: Number(e.target.value) })
            }
          />
        </div>
        <div className="control-group">
          <label htmlFor="resolution">
            Resolution: {docState.network.resolution}
          </label>
          <input
            id="resolution"
            type="range"
            min="0.1"
            max="2"
            step="0.1"
            value={docState.network.resolution}
            onChange={(e) =>
              updateSettings({ resolution: Number(e.target.value) })
            }
          />
        </div>
        {selectedNode && (
          <div className="node-info">
            <h4>Selected Node</h4>
            <p>ID: {selectedNode}</p>
            {networkData && (
              <>
                <p>
                  Cluster:{" "}
                  {
                    networkData.nodes.find((n) => n.id === selectedNode)
                      ?.cluster
                  }
                </p>
                <p>
                  Degree:{" "}
                  {networkData.nodes.find((n) => n.id === selectedNode)?.degree}
                </p>
              </>
            )}
          </div>
        )}
      </div>
    </>
  );
};
