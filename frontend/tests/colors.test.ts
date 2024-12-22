import { describe, expect, test } from "bun:test";
import { type ColorScaleArray, colorScales } from "../src/colorScales";
import {
  findScaleLower,
  findScaleUpper,
  interpolateColor,
  originalRgbFormat,
} from "../src/colors";

const testScale: ColorScaleArray = [
  [0, "rgb(8,29,88)"],
  [0.125, "rgb(37,52,148)"],
  [0.25, "rgb(34,94,168)"],
  [0.375, "rgb(29,145,192)"],
  [0.5, "rgb(65,182,196)"],
  [0.625, "rgb(127,205,187)"],
  [0.75, "rgb(199,233,180)"],
  [0.875, "rgb(237,248,217)"],
  [1, "rgb(255,255,217)"],
];

test("findScaleLower works", () => {
  expect(findScaleLower(testScale, 0.25)).toEqual([0.25, "rgb(34,94,168)"]);
});
test("findScaleUpper works", () => {
  expect(findScaleUpper(testScale, 0.9)).toEqual([1, "rgb(255,255,217)"]);
});

describe("interpolateColor", () => {
  test("works with a small dataset", () => {
    const interpolatedColor = interpolateColor(
      [
        [0, "rgb(0, 0, 0)"],
        [1, "rgb(255, 255, 255)"],
      ],
      0.5,
      originalRgbFormat,
    );
    expect(interpolatedColor).toEqual({
      value: [0.5, "rgb(128, 128, 128)"],
      lower: [0, "rgb(0, 0, 0)"],
      upper: [1, "rgb(255, 255, 255)"],
    });
  });

  describe("with a real dataset", () => {
    test("works with normal ratio", () => {
      const interpolatedColor = interpolateColor(
        colorScales.Yellow_Green_Blue,
        0.57,
        originalRgbFormat,
      );
      expect(interpolatedColor).toEqual({
        value: [0.5599999999999996, "rgb(100, 195, 191)"],
        lower: [0.5, "rgb(65,182,196)"],
        upper: [0.625, "rgb(127,205,187)"],
      });
    });

    test("works with zero ratio", () => {
      const interpolatedColor = interpolateColor(
        colorScales.Yellow_Green_Blue,
        0.5,
        originalRgbFormat,
      );
      expect(interpolatedColor).toEqual({
        value: [0, "rgb(65,182,196)"],
        lower: [0.5, "rgb(65,182,196)"],
        upper: [0.5, "rgb(65,182,196)"],
      });
    });
  });
});
