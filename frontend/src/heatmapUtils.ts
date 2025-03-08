export const getFontSizeForCell = (
  cellWidth: number,
  characterCount: number,
) => {
  const characterWidth = 0.6; // Approximate width for Roboto Mono
  const availableWidth = cellWidth * 0.8;

  return Math.min(
    20,
    Math.max(
      0.25,
      Math.floor(availableWidth / (characterCount * characterWidth)),
    ),
  );
};
