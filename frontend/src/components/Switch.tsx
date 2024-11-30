import type React from "react";
import { Switch as RACSwitch, type SwitchProps } from "react-aria-components";

export const Switch = ({
  children,
  ...props
}: React.PropsWithChildren<SwitchProps>) => {
  return (
    <>
      <RACSwitch {...props}>
        <div className="indicator" />
        {children}
      </RACSwitch>
    </>
  );
};
