import {
  OverlayArrow,
  Tooltip as RACTooltip,
  TooltipTrigger,
  type TooltipTriggerComponentProps,
} from "react-aria-components";

export const Tooltip = ({
  tooltip,
  children,
  ...props
}: React.PropsWithChildren<
  TooltipTriggerComponentProps & { tooltip: React.ReactNode }
>) => (
  <TooltipTrigger {...props}>
    {children}
    <RACTooltip>
      <OverlayArrow>
        <svg aria-hidden="true" width={8} height={8} viewBox="0 0 8 8">
          <path d="M0 0 L4 4 L8 0" />
        </svg>
      </OverlayArrow>
      {tooltip}
    </RACTooltip>
  </TooltipTrigger>
);
