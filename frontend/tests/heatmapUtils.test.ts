import { describe, expect, it } from "bun:test";
import type { ColorScaleArray } from "../src/colorScales";
import {
  formatClustermapData,
  formatHeatmapData,
  getCellMetrics,
} from "../src/heatmapUtils";
import type { HeatmapSettings } from "../src/plotTypes";

describe("getCellMetrics", () => {
  it("should return the correct metrics", () => {
    expect(getCellMetrics(20, 1, 5)).toEqual({
      cellSize: 19.5,
      cellSpace: 0.5,
      cellOffset: 0.25,
      fontSize: 5.2,
      textOffset: 9.75,
    });
  });
  it("should throw an error if characterCount is less than or equal to 0", () => {
    expect(() => getCellMetrics(10, 2, 0)).toThrowError(
      "characterCount must be > 0",
    );
    expect(() => getCellMetrics(10, 2, -1)).toThrowError(
      "characterCount must be > 0",
    );
  });
});

describe("formatHeatmapData", () => {
  it("should return formatted data", () => {
    const data = [
      [100, null, null, null],
      [72.6, 100, null, null],
      [72.62, 97.61, 100, null],
      [8.2, 66.99, 66.88, 100],
    ];

    const settings = {
      colorScaleKey: "Test",
      vmax: 100,
      vmin: 8.2,
      annotation_rounding: 2,
    };

    const colorScale: ColorScaleArray = [
      [0, "rgb(0, 0, 0)"],
      [1, "rgb(255, 255, 255)"],
    ];

    expect(
      formatHeatmapData(data, settings as HeatmapSettings, colorScale),
    ).toMatchSnapshot();
  });
});

describe("formatClustermapData", () => {
  it("should return formatted data", () => {
    const data = [
      [100, null, null, null],
      [72.6, 100, null, null],
      [72.62, 97.61, 100, null],
      [8.2, 66.99, 66.88, 100],
    ];

    const tickText = ["A", "B", "C", "D"];

    const clusterData = [
      { id: "A", group: 1 },
      { id: "B", group: 2 },
      { id: "C", group: 2 },
      { id: "D", group: 3 },
    ];

    expect(formatClustermapData(data, tickText, clusterData)).toMatchSnapshot();
  });
});
