import React from 'react';
import { createRoot } from 'react-dom/client';
import './index.scss';
import { App } from './components/App';

const container = document.getElementById('app');
const root = createRoot(container);
root.render(<App />);
