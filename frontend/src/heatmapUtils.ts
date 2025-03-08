export const getPlotMetrics = (width: number, height: number) => {
  // WIP
  return {
    width: width,
    height: height,
  };
};

export const getCellMetrics = (
  cellSize: number,
  cellSpace: number,
  characterCount: number,
) => {
  if (characterCount <= 0) throw new Error("characterCount must be > 0");

  const CHARACTER_WIDTH = 0.6; // Approximate width for Roboto Mono @ 10px
  const USABLE_SPACE_RATIO = 0.8;
  const MIN_FONT_SIZE = 0.25;
  const MAX_FONT_SIZE = 20;

  const scaledCellSpace = cellSpace * (cellSize / (cellSize + 20));
  const spacedCellSize = Math.max(1, cellSize - scaledCellSpace);
  const scaledFontSize = Math.min(
    MAX_FONT_SIZE,
    Math.max(
      MIN_FONT_SIZE,
      (spacedCellSize * USABLE_SPACE_RATIO) /
        (characterCount * CHARACTER_WIDTH),
    ),
  );

  return {
    cellSize: spacedCellSize,
    cellSpace: scaledCellSpace,
    cellOffset: scaledCellSpace / 2,
    fontSize: scaledFontSize,
    textOffset: spacedCellSize / 2,
  };
};
