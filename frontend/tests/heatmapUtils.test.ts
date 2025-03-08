import { describe, expect, it } from "bun:test";
import { getFontSizeForCell } from "../src/heatmapUtils";

describe("getFontSizeForCell", () => {
  it("should return the correct font size for a cell", () => {
    expect(getFontSizeForCell(100, 5)).toBe(26);
    expect(getFontSizeForCell(50, 5)).toBe(13);
    expect(getFontSizeForCell(25, 3)).toBe(11);
    expect(getFontSizeForCell(1, 3)).toBe(1);
    expect(getFontSizeForCell(0, 1)).toBe(1);
  });
});
