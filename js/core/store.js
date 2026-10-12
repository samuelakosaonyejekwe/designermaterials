// Persistence (IndexedDB with in-memory fallback) and the project model: shared data bus + run records.
const mem = new Map(); let dbp = null;
function db() {
  if (dbp) return dbp;
  dbp = new Promise(res => { try { const r = indexedDB.open('designer-materials', 1); r.onupgradeneeded = () => r.result.createObjectStore('kv'); r.onsuccess = () => res(r.result); r.onerror = () => res(null); } catch { res(null); } });
  return dbp;
}
const tx = async (mode, fn) => { const d = await db(); if (!d) return fn(null); return new Promise((res, rej) => { const t = d.transaction('kv', mode), rq = fn(t.objectStore('kv')); t.oncomplete = () => res(rq && rq.result); t.onerror = () => rej(t.error); }); };
export const kv = {
  get: async k => { try { const v = await tx('readonly', s => s ? s.get(k) : null); return v === undefined || v === null ? mem.get(k) : v; } catch { return mem.get(k); } },
  set: async (k, v) => { mem.set(k, v); try { await tx('readwrite', s => s && s.put(v, k)); } catch { /* quota or private mode: memory only */ } },
  del: async k => { mem.delete(k); try { await tx('readwrite', s => s && s.delete(k)); } catch { /* ignore */ } },
  keys: async () => { try { return (await tx('readonly', s => s ? s.getAllKeys() : null)) || [...mem.keys()]; } catch { return [...mem.keys()]; } },
};
export const prefs = {
  get(k, d) { try { const v = localStorage.getItem('dm.' + k); return v === null ? d : JSON.parse(v); } catch { return d; } },
  set(k, v) { try { localStorage.setItem('dm.' + k, JSON.stringify(v)); } catch { /* ignore */ } },
};

// ---- project ----
export const project = { id: 'default', name: 'Untitled project', bus: {}, runs: {}, inputs: {}, geometry: [], datasets: [], notes: '' };
const listeners = new Set();
export const onProject = fn => { listeners.add(fn); return () => listeners.delete(fn); };
const emit = what => listeners.forEach(f => { try { f(what); } catch (e) { console.error(e); } });
let saveT = null;
export function saveProject() { clearTimeout(saveT); saveT = setTimeout(() => kv.set('project:' + project.id, JSON.parse(JSON.stringify(project, (k, v) => ArrayBuffer.isView(v) ? Array.from(v) : v))), 250); }
export async function loadProject(id) {
  id = id || prefs.get('project', 'default'); const p = await kv.get('project:' + id);
  Object.assign(project, { id, name: 'Untitled project', bus: {}, runs: {}, inputs: {}, geometry: [], datasets: [], notes: '' }, p || {}); prefs.set('project', id); emit('load'); return project;
}
export async function listProjects() { const out = []; for (const k of await kv.keys()) if (String(k).startsWith('project:')) { const p = await kv.get(k); if (p) out.push({ id: p.id, name: p.name, runs: Object.keys(p.runs || {}).length }); } return out; }
export async function newProject(name) { const id = 'p' + Date.now().toString(36); Object.assign(project, { id, name: name || 'New project', bus: {}, runs: {}, inputs: {}, geometry: [], datasets: [], notes: '' }); prefs.set('project', id); saveProject(); emit('load'); return project; }
export async function deleteProject(id) { await kv.del('project:' + id); }
export function busValues() { const o = {}; for (const k in project.bus) o[k] = project.bus[k].value; return o; }
export function publish(suiteId, values, when) { for (const k in values || {}) project.bus[k] = { value: values[k], from: suiteId, at: when }; saveProject(); emit('bus'); }
export function setBus(key, value, from = 0) { project.bus[key] = { value, from, at: new Date().toISOString() }; saveProject(); emit('bus'); }
export function recordRun(suiteId, rec) { project.runs[suiteId] = rec; saveProject(); emit('run'); }
export function setInputs(suiteId, values) { project.inputs[suiteId] = values; saveProject(); }
export function exportProject() { return JSON.stringify({ format: 'designer-materials-project', version: 1, exported: new Date().toISOString(), project }, (k, v) => ArrayBuffer.isView(v) ? Array.from(v) : v, 1); }
export function importProject(text) { const o = JSON.parse(text); if (o.format !== 'designer-materials-project') throw new Error('Not a Designer Materials project file'); Object.assign(project, o.project, { id: 'p' + Date.now().toString(36) }); prefs.set('project', project.id); saveProject(); emit('load'); }
