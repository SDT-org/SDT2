import React from "react";
import tinycolor from "tinycolor2";
import type { AppState } from "../appState";
import type { ColorScaleArray } from "../colorScales";
import { interpolateColor, originalRgbFormat } from "../colors";
import type { HeatmapData } from "../plotTypes";

export const useAnnotations = (
  enabled: boolean,
  data: HeatmapData,
  vmin: number,
  vmax: number,
  _colorScale: ColorScaleArray,
  reverse: boolean,
  roundTo: AppState["client"]["heatmap"]["annotation_rounding"],
) =>
  React.useMemo(() => {
    if (!enabled) {
      return {
        x: [],
        y: [],
        text: [],
        textColors: [],
      };
    }

    const x: number[] = [];
    const y: number[] = [];
    const text: string[] = [];
    const textColors: (string | null)[] = [];
    let colorScale = _colorScale;

    if (reverse) {
      colorScale = [..._colorScale]
        .reverse()
        .map((data, i) => [
          (_colorScale[i] ?? _colorScale[0])[0],
          data[1],
        ]) as ColorScaleArray;
    }

    const dataMin = vmin;
    const dataMax = vmax;
    const dataDiff = dataMax - dataMin;

    data.forEach((row, rowIndex) => {
      row.forEach((datum, columnIndex) => {
        x.push(columnIndex);
        y.push(rowIndex);

        if (datum === null) {
          text.push("");
        } else {
          text.push(Number.parseFloat(datum).toFixed(roundTo));
        }

        const parsedDatum = Number.parseFloat(datum);
        const normalizedDatum = Math.max(0, (parsedDatum - dataMin) / dataDiff);

        if (datum === null) {
          textColors.push(null);
        } else {
          const interpolated = interpolateColor(
            colorScale,
            normalizedDatum,
            originalRgbFormat,
          );
          const textColor = tinycolor
            .mostReadable(interpolated.value[1], ["#fff", "#000"], {
              includeFallbackColors: false,
              level: "AAA",
              size: "small",
            })
            .toHexString();

          textColors.push(textColor);
        }
      });
    });

    return {
      x,
      y,
      text,
      textColors,
    };
  }, [enabled, data, vmin, vmax, _colorScale, reverse, roundTo]);
