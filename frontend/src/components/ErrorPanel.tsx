import { TabPanel } from "react-aria-components";
import type { DocState } from "../appState";

export const ErrorPanel = ({ docState }: { docState: DocState }) => (
  <TabPanel
    id={docState.id}
    key={docState.id}
    className="app-panel app-panel-full-width"
  >
    <div className="app-main centered error">
      <h2>An error occured while processing this file.</h2>
      <p>Please ensure it is a valid, SDT-compatible file.</p>
      <details>
        <summary>Error details</summary>
        <pre>{docState.invalid?.reason}</pre>
      </details>
    </div>
  </TabPanel>
);
