import type {
  ColorPickerProps as RACColorPickerProps,
  ColorSliderProps as RACColorSliderProps,
} from "react-aria-components";
import {
  Button,
  ColorArea,
  ColorSwatch,
  ColorThumb,
  Dialog,
  DialogTrigger,
  Label,
  Popover,
  ColorPicker as RACColorPicker,
  ColorSlider as RACColorSlider,
  SliderOutput,
  SliderTrack,
} from "react-aria-components";

interface ColorSliderProps extends RACColorSliderProps {
  label?: string;
}

export const ColorSlider = ({ label, ...props }: ColorSliderProps) => {
  return (
    <RACColorSlider {...props}>
      <Label>{label}</Label>
      <SliderOutput />
      <SliderTrack
        style={({ defaultStyle }) => ({
          background: `${defaultStyle.background},
            repeating-conic-gradient(#CCC 0% 25%, white 0% 50%) 50% / 16px 16px`,
        })}
      >
        <ColorThumb />
      </SliderTrack>
    </RACColorSlider>
  );
};

interface ColorPickerProps extends RACColorPickerProps {
  label?: string;
  children?: React.ReactNode;
  value: string;
}

export const ColorPicker = ({
  label,
  children,
  value,
  ...props
}: ColorPickerProps) => {
  return (
    <RACColorPicker {...props}>
      <DialogTrigger>
        <Button className="color-picker">
          <ColorSwatch color={value} />
          <span>{label}</span>
        </Button>
        <Popover placement="bottom start">
          <Dialog className="color-picker-dialog">
            {children || (
              <>
                <ColorArea
                  colorSpace="hsb"
                  xChannel="saturation"
                  yChannel="brightness"
                  value={value}
                >
                  <ColorThumb />
                </ColorArea>
                <ColorSlider colorSpace="hsl" channel="hue" value={value} />
                <ColorSlider colorSpace="hsl" channel="alpha" value={value} />
              </>
            )}
          </Dialog>
        </Popover>
      </DialogTrigger>
    </RACColorPicker>
  );
};
