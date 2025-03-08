import { describe, expect, it } from "bun:test";
import { getCellMetrics } from "../src/heatmapUtils";

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
