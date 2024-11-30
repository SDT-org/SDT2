import { createRoot } from "react-dom/client";
import "./index.scss";
import { App } from "./components/App";

const container = document.getElementById("app");
if (!container) {
  throw new Error("could not find container to mount app");
}

const root = createRoot(container);
root.render(<App />);
