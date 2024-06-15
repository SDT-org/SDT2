import React from 'react';
import { AppState } from '../appState';

export const Loader = ({
  appState,
}: {
  appState: AppState;
}) =>
  <div className="app-wrapper">
    <div className="app-main centered">
      <h4>Analyzing...</h4>
      <progress id="progress" max="100" value={appState.progress}>
        {appState.progress} %
      </progress>
    </div>
  </div>
