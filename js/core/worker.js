// Background solver: runs a suite off the main thread so long, full-resolution simulations keep the interface responsive.
let cancelled = false;
self.onmessage = async e => {
  const m = e.data;
  if (m.type === 'cancel') { cancelled = true; return; }
  if (m.type !== 'run') return;
  cancelled = false;
  try {
    const suite = (await import(`../suites/s${String(m.id).padStart(2, '0')}.js`)).default, t0 = performance.now();
    const result = await suite.run(m.inputs, { bus: m.bus || {}, hifi: m.hifi, cancelled: () => cancelled, onProgress: p => self.postMessage({ type: 'progress', ...p }) });
    self.postMessage({ type: 'done', result, runtime: performance.now() - t0 });
  } catch (err) { self.postMessage({ type: 'error', message: String(err && err.message || err) }); }
};
