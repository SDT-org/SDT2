export const getScaledFontSize = (base: number, count: number) =>
  Math.min(16, Math.max(1, base / (1 + count / 60)));
