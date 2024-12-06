import React from "react";

export const useAppBlur = () =>
  React.useEffect(() => {
    const addBlur = () => {
      document.body.classList.add("blur");
    };

    const removeBlur = () => {
      document.body.classList.remove("blur");
    };

    window.addEventListener("focus", removeBlur);
    window.addEventListener("blur", addBlur);

    return () => {
      window.removeEventListener("focus", removeBlur);
      window.removeEventListener("blur", addBlur);
    };
  }, []);
