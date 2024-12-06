import React from "react";

export const useWaitForPywebview = (callback: () => void) =>
  React.useEffect(() => {
    const waitForPywebview = () =>
      new Promise((resolve) => {
        if (!window.pywebview) {
          setTimeout(() => resolve(waitForPywebview()), 10);
        } else {
          resolve(true);
        }
      });

    waitForPywebview().then(callback);
  }, [callback]);
