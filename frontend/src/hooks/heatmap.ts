import React from "react";
import type { ClustermapSettings, HeatmapSettings } from "../plotTypes";

export const useSize = (
  elementRef: React.MutableRefObject<HTMLDivElement | null>,
  leftSidebarCollapsed: boolean,
) => {
  const [size, setSize] = React.useState<{ width: number; height: number }>({
    width: 0,
    height: 0,
  });

  const updateSize = React.useCallback(() => {
    if (elementRef.current) {
      const { offsetWidth, offsetHeight } = elementRef.current;
      setSize({ width: offsetWidth, height: offsetHeight });
    }
  }, [elementRef.current]);

  // TODO: Do this the right way
  // biome-ignore lint/correctness/useExhaustiveDependencies(leftSidebarCollapsed): trigger updateSize
  React.useEffect(() => {
    updateSize();

    const handleResize = () => {
      updateSize();
    };

    window.addEventListener("resize", handleResize);

    return () => {
      window.removeEventListener("resize", handleResize);
    };
  }, [updateSize, leftSidebarCollapsed]);

  return size;
};

export const useMetrics = (
  settings: HeatmapSettings | ClustermapSettings,
  tickText: string[],
) => {
  const longestTickWidth = React.useMemo(
    () =>
      Math.max(...tickText.map((tick) => tick.length)) *
      settings.axlabel_fontsize,
    [tickText, settings.axlabel_fontsize],
  );

  const margin = React.useMemo(
    () => ({
      top: 60,
      right: 60,
      bottom: settings.axis_labels ? Math.max(longestTickWidth, 60) : 60,
      left: settings.axis_labels ? Math.max(longestTickWidth, 60) : 60,
    }),
    [longestTickWidth, settings.axis_labels],
  );

  return {
    ...("cbar_shrink" in settings && {
      cbar_shrink: settings.cbar_shrink * 60,
    }),
    ...("cbar_aspect" in settings && {
      cbar_aspect: settings.cbar_aspect * 10,
    }),
    margin,
  };
};
