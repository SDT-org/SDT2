import { defineConfig } from "vite";
import react from "@vitejs/plugin-react";

export default defineConfig({
  plugins: [
    react(),
    commonjs({
      include: ["frontend/vendor/plotly-custom.min.js"],
    }),
  ],
  build: {
    outDir: "gui",
    chunkSizeWarningLimit: 5000,
  },
});
