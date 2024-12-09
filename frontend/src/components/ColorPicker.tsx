import type {
  ColorPickerProps as RACColorPickerProps,
  ColorSliderProps as RACColorSliderProps,
} from "react-aria-components";
import {
  Button,
  ColorArea,
  ColorSwatch,
  ColorSwatchPicker,
  ColorSwatchPickerItem,
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
import { Colors } from "../colors";

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
  className?: string;
  "aria-label"?: string;
}

export const ColorPicker = ({
  label,
  children,
  value,
  className,
  ...props
}: ColorPickerProps) => {
  return (
    <RACColorPicker {...props}>
      <DialogTrigger>
        <Button
          className={`color-picker ${className || ""}`.trimEnd()}
          {...(props["aria-label"]
            ? { "aria-label": props["aria-label"] }
            : {})}
        >
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
                <ColorSwatchPicker>
                  <ColorSwatchPickerItem color="#3a5975">
                    <ColorSwatch />
                  </ColorSwatchPickerItem>
                  <ColorSwatchPickerItem color="#467588">
                    <ColorSwatch />
                  </ColorSwatchPickerItem>
                  <ColorSwatchPickerItem color="#519398">
                    <ColorSwatch />
                  </ColorSwatchPickerItem>
                  <ColorSwatchPickerItem color="#5bb3a1">
                    <ColorSwatch />
                  </ColorSwatchPickerItem>
                  <ColorSwatchPickerItem color="#62d49e">
                    <ColorSwatch />
                  </ColorSwatchPickerItem>
                </ColorSwatchPicker>
                <ColorSwatchPicker>
                  <ColorSwatchPickerItem color="#e8960a">
                    <ColorSwatch />
                  </ColorSwatchPickerItem>
                  <ColorSwatchPickerItem color="#d18813">
                    <ColorSwatch />
                  </ColorSwatchPickerItem>
                  <ColorSwatchPickerItem color="#df6056">
                    <ColorSwatch />
                  </ColorSwatchPickerItem>
                  <ColorSwatchPickerItem color="#da3f42">
                    <ColorSwatch />
                  </ColorSwatchPickerItem>
                  <ColorSwatchPickerItem color="#d20c31">
                    <ColorSwatch />
                  </ColorSwatchPickerItem>
                </ColorSwatchPicker>
              </>
            )}
          </Dialog>
        </Popover>
      </DialogTrigger>
    </RACColorPicker>
  );
};
