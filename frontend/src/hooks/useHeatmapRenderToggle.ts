import React from "react";

export const useHeatmapRenderToggle = () => {
  const [forceSvgRender, setForceSvgRender] = React.useState(false);

  React.useEffect(() => {
    const handleKeyDown = (event: KeyboardEvent) => {
      if ((event.metaKey || event.altKey) && event.key === "1") {
        setForceSvgRender(true);
        event.preventDefault();
      } else if ((event.metaKey || event.altKey) && event.key === "2") {
        setForceSvgRender(false);
        event.preventDefault();
      }
    };

    document.addEventListener("keydown", handleKeyDown);
    return () => document.removeEventListener("keydown", handleKeyDown);
  }, []);

  return forceSvgRender;
};
