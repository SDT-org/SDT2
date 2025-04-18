export const getScaledFontSize = (base: number, count: number) =>
  count < 10 ? 16 : Math.min(16, Math.max(1, base / (1 + count / 60)));
