import React from "react";
import {
  Input,
  Label,
  LabelContext,
  NumberField,
  Slider as RACSlider,
  type SliderProps as RACSliderProps,
  SliderOutput,
  SliderStateContext,
  SliderThumb,
  SliderTrack,
  useSlottedContext,
} from "react-aria-components";

interface SliderProps<T> extends RACSliderProps<T> {
  label?: string;
  labelClassName?: string;
  thumbLabels?: string[];
  includeField?: boolean;
}

const SliderNumberField = () => {
  const state = React.useContext(SliderStateContext);
  const labelProps = useSlottedContext(LabelContext);
  return (
    <NumberField
      aria-labelledby={labelProps?.id ? labelProps.id : ""}
      value={state?.values?.[0] ? state.values[0] : 0}
      onChange={(v) => (state ? state.setThumbValue(0, v) : null)}
    >
      <Input />
    </NumberField>
  );
};

export const Slider = <T extends number | number[]>({
  label,
  thumbLabels,
  labelClassName,
  includeField,
  ...props
}: SliderProps<T>) => {
  return (
    <RACSlider {...props}>
      {label && <Label className={labelClassName}>{label}</Label>}
      {includeField ? (
        <SliderNumberField />
      ) : (
        <SliderOutput>
          {({ state }) => <>{state.getThumbValueLabel(0)}</>}
        </SliderOutput>
      )}
      <SliderTrack>
        {({ state }) => (
          <>
            <div className="track" />
            <div
              className="fill"
              style={{
                width: `${state.getThumbPercent(0) * 100}%`,
              }}
            />
            <SliderThumb />
          </>
        )}
      </SliderTrack>
    </RACSlider>
  );
};
