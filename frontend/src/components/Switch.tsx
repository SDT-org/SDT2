import React from "react";
import { Switch as RACSwitch, SwitchProps } from "react-aria-components";

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
