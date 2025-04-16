# TypeScript Fixes for Plot Components

## Changes Made

1. Removed "Plot Titles" references in all sidebar components and changed to "Plot Title"
2. Removed X and Y axis title fields from ViolinSidebar and HistogramSidebar components
3. Removed X and Y axis title rendering code from Violin and Histogram components
4. Fixed TypeScript errors:

### Histogram.tsx
- Changed `tickfont: { weight: "normal" }` to `tickfont: { family: "sans-serif", size: 11 }`
- Removed the `@ts-ignore` comment and replaced `weight: "bold"` with `size: 14`

### Violin.tsx
- Fixed tickmode typing by adding `as const` to the string literal type: `tickmode: "array" as const`
- Replaced `...plotFontMonospace` with explicit property assignment: `family: plotFontMonospace.family, size: 11`
- Moved `boxgap` property into the xaxis configuration (was incorrectly at top level)
- Removed `weight` property from font objects

### useRelayoutUpdateTitles.ts
- Removed axis title update functionality since we no longer need it
- Updated type definitions to remove axis title properties

### distributionState.ts
- Removed `xtitle` and `ytitle` properties from the violin and histogram visualization types