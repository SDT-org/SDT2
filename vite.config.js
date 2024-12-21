import { defineConfig } from "vite";
import react from "@vitejs/plugin-react";

export default defineConfig(({ mode }) => ({
  root: 'frontend/src',
  plugins: [
    react(),
    mode === 'development' && {
     name: 'inject-devtools-script',
     transformIndexHtml(html) {
       return html.replace(
         '</body>',
         '<script src="http://localhost:8097"></script></body>'
       );
     },
   },
  ],
  build: {
    outDir: "../../gui",
    chunkSizeWarningLimit: 5000,
    emptyOutDir: true,
  },
  publicDir: "../../docs",
}));
