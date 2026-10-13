// Entry point of the single-file edition: the same code runs as the page and as its background worker.
if (typeof document === 'undefined') await import('./core/worker.js'); else await import('./app.js');
