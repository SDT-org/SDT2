import { defineConfig } from "vite";
import react from "@vitejs/plugin-react";

export default defineConfig({
  root: 'frontend/src',
  plugins: [
    react(),
  ],
  build: {
    outDir: "../../gui",
    chunkSizeWarningLimit: 5000,
    emptyOutDir: true,
  },
  publicDir: "../../docs"
});
