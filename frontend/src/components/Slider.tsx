import {
  Label,
  Slider as RACSlider,
  type SliderProps as RACSliderProps,
  SliderOutput,
  SliderThumb,
  SliderTrack,
} from "react-aria-components";

import React from "react";

interface SliderProps<T> extends RACSliderProps<T> {
  label?: string;
  labelClassName?: string;
  thumbLabels?: string[];
}

export const Slider = <T extends number | number[]>({
  label,
  thumbLabels,
  labelClassName,
  ...props
}: SliderProps<T>) => {
  // TODO: .fill Doesn't work with multiple values yet
  return (
    <RACSlider {...props}>
      {label && <Label className={labelClassName}>{label}</Label>}
      <SliderOutput>
        {({ state }) =>
          state.values.map((_, i) => state.getThumbValueLabel(i)).join(" – ")
        }
      </SliderOutput>
      <SliderTrack>
        {({ state }) => (
          <>
            <div className="track" />
            <div
              className="fill"
              style={{
                width: state.getThumbPercent(0) * 100 + "%",
              }}
            />
            {state.values.map((_, i) => (
              <SliderThumb
                key={i}
                index={i}
                aria-label={thumbLabels?.[i] ?? ""}
              />
            ))}
          </>
        )}
      </SliderTrack>
    </RACSlider>
  );
};
