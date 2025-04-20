import React from "react";

export const useWaitForPywebview = (callback: () => void) => {
  const hasRunRef = React.useRef(false);

  React.useEffect(() => {
    let checkInterval: number | undefined = undefined;

    if (hasRunRef.current) {
      return;
    }

    const checkForPywebview = () => {
      if (window.pywebview) {
        clearInterval(checkInterval);
        hasRunRef.current = true;
        callback();
      }
    };

    checkInterval = window.setInterval(checkForPywebview, 10);

    return () => {
      if (checkInterval) {
        clearInterval(checkInterval);
      }
    };
  }, [callback]);
};
