import { align } from './core.js';
self.onmessage = ({ data }) => {
  try {
    const started = performance.now(),
      result = align(data.reference, data.query, data.options);
    self.postMessage({ result, milliseconds: performance.now() - started }, [result.matrix.buffer]);
  } catch (error) {
    self.postMessage({ error: error.message });
  }
};
