import {
  Label,
  Slider as RACSlider,
  type SliderProps as RACSliderProps,
  SliderOutput,
  SliderThumb,
  SliderTrack,
} from "react-aria-components";

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
  return (
    <RACSlider {...props}>
      {label && <Label className={labelClassName}>{label}</Label>}
      <SliderOutput>
        {({ state }) => <>{state.getThumbValueLabel(0)}</>}
      </SliderOutput>
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
