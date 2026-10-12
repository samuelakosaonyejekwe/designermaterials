// Chemistry kernel: periodic table, formula parsing, molar mass, equation balancing (mass + charge),
// molecular / complete-ionic / net-ionic equation generation and pretty-printing.

// symbol:mass:electronegativity:covalentRadius(Å)
const RAW = 'H:1.008:2.20:0.31 He:4.0026:0:0.28 Li:6.94:0.98:1.28 Be:9.0122:1.57:0.96 B:10.81:2.04:0.84 C:12.011:2.55:0.76 N:14.007:3.04:0.71 O:15.999:3.44:0.66 F:18.998:3.98:0.57 Ne:20.180:0:0.58 Na:22.990:0.93:1.66 Mg:24.305:1.31:1.41 Al:26.982:1.61:1.21 Si:28.085:1.90:1.11 P:30.974:2.19:1.07 S:32.06:2.58:1.05 Cl:35.45:3.16:1.02 Ar:39.948:0:1.06 K:39.098:0.82:2.03 Ca:40.078:1.00:1.76 Sc:44.956:1.36:1.70 Ti:47.867:1.54:1.60 V:50.942:1.63:1.53 Cr:51.996:1.66:1.39 Mn:54.938:1.55:1.39 Fe:55.845:1.83:1.32 Co:58.933:1.88:1.26 Ni:58.693:1.91:1.24 Cu:63.546:1.90:1.32 Zn:65.38:1.65:1.22 Ga:69.723:1.81:1.22 Ge:72.630:2.01:1.20 As:74.922:2.18:1.19 Se:78.971:2.55:1.20 Br:79.904:2.96:1.20 Kr:83.798:3.00:1.16 Rb:85.468:0.82:2.20 Sr:87.62:0.95:1.95 Y:88.906:1.22:1.90 Zr:91.224:1.33:1.75 Nb:92.906:1.60:1.64 Mo:95.95:2.16:1.54 Tc:98:1.90:1.47 Ru:101.07:2.20:1.46 Rh:102.91:2.28:1.42 Pd:106.42:2.20:1.39 Ag:107.87:1.93:1.45 Cd:112.41:1.69:1.44 In:114.82:1.78:1.42 Sn:118.71:1.96:1.39 Sb:121.76:2.05:1.39 Te:127.60:2.10:1.38 I:126.90:2.66:1.39 Xe:131.29:2.60:1.40 Cs:132.91:0.79:2.44 Ba:137.33:0.89:2.15 La:138.91:1.10:2.07 Ce:140.12:1.12:2.04 Pr:140.91:1.13:2.03 Nd:144.24:1.14:2.01 Pm:145:1.13:1.99 Sm:150.36:1.17:1.98 Eu:151.96:1.20:1.98 Gd:157.25:1.20:1.96 Tb:158.93:1.10:1.94 Dy:162.50:1.22:1.92 Ho:164.93:1.23:1.92 Er:167.26:1.24:1.89 Tm:168.93:1.25:1.90 Yb:173.05:1.10:1.87 Lu:174.97:1.27:1.87 Hf:178.49:1.30:1.75 Ta:180.95:1.50:1.70 W:183.84:2.36:1.62 Re:186.21:1.90:1.51 Os:190.23:2.20:1.44 Ir:192.22:2.20:1.41 Pt:195.08:2.28:1.36 Au:196.97:2.54:1.36 Hg:200.59:2.00:1.32 Tl:204.38:1.62:1.45 Pb:207.2:2.33:1.46 Bi:208.98:2.02:1.48 Po:209:2.00:1.40 At:210:2.20:1.50 Rn:222:2.20:1.50 Fr:223:0.79:2.60 Ra:226:0.90:2.21 Ac:227:1.10:2.15 Th:232.04:1.30:2.06 Pa:231.04:1.50:2.00 U:238.03:1.38:1.96 Np:237:1.36:1.90 Pu:244:1.28:1.87 Am:243:1.30:1.80 Cm:247:1.30:1.69';
export const ELEMENTS = {}; export const SYMBOLS = [];
RAW.split(' ').forEach((t, i) => { const [s, m, en, r] = t.split(':'); ELEMENTS[s] = { Z: i + 1, symbol: s, mass: +m, en: +en, radius: +r }; SYMBOLS.push(s); });
const CPK = { H: '#e8e8e8', C: '#555', N: '#3050f8', O: '#ff2d2d', F: '#90e050', S: '#e6c200', P: '#ff8000', Cl: '#1ff01f', Zn: '#7d80b0', Ti: '#bfc2c7', Si: '#f0c8a0', Au: '#ffd123', Ag: '#c0c0c0', Cu: '#c88033', Fe: '#e06633', Al: '#bfa6a6', Ga: '#c28f8f', As: '#bd80e3', Cd: '#ffd98f', Se: '#ffa100', Mo: '#54b5b5', W: '#2194d6', Zr: '#94e0e0', Ba: '#00c900', Sr: '#00ff00', Pb: '#575961', Ni: '#50d050', Co: '#f090a0', Mn: '#9c7ac7', Li: '#cc80ff', Na: '#ab5cf2', K: '#8f40d4', Mg: '#8aff00', Ca: '#3dff00', B: '#ffb5b5', In: '#a67573', Sn: '#668080', Pt: '#d0d0e0', Ce: '#ffffc7', V: '#a6a6ab', Nb: '#73c2c9', Te: '#d47a00', Bi: '#9e4fb5', Ge: '#668f8f', Y: '#94ffff', Hf: '#4dc2ff', Sc: '#e6e6e6' };
export const elementColor = s => CPK[s] || '#b0a0c0';

const SUB = '₀₁₂₃₄₅₆₇₈₉', SUPD = '⁰¹²³⁴⁵⁶⁷⁸⁹';
export const subscript = s => String(s).replace(/\d/g, d => SUB[d]).replace('.', '.');
export const superscript = s => String(s).replace(/\d/g, d => SUPD[d]).replace('+', '⁺').replace('-', '⁻');

// Charge convention: "SO4^2-" (caret, explicit); "Zn2+", "Fe3+", "S2-" (single element: digits are the charge);
// "[Zn(OH)4]2-" (bracketed complex: digits are the charge); "NO3-", "NH4+", "OH-" (polyatomic: digits are stoichiometric, charge ±1).
function splitCharge(f) {
  let m = f.match(/\^(\d*)([+-])$/);
  if (m) return { body: f.slice(0, m.index), charge: (m[1] ? +m[1] : 1) * (m[2] === '+' ? 1 : -1) };
  m = f.match(/^(.*?)(\d*)([+-])$/);
  if (!m || !m[1]) return { body: f, charge: 0 };
  const sign = m[3] === '+' ? 1 : -1;
  if (/^[A-Z][a-z]?$/.test(m[1]) || /[\]}]$/.test(m[1])) return { body: m[1], charge: (m[2] ? +m[2] : 1) * sign };
  return { body: m[1] + m[2], charge: sign };
}
// Parse "Zn(NO3)2·6H2O", "SO4^2-", "Zn2+", "[Zn(OH)4]2-", "e-". Returns {atoms:{El:n}, charge, phase}
export function parseFormula(input) {
  let f = String(input).trim(); let phase = '';
  const ph = f.match(/\((s|l|g|aq)\)$/); if (ph) { phase = ph[1]; f = f.slice(0, ph.index); }
  f = f.replace(/[₀-₉]/g, c => SUB.indexOf(c)).replace(/[⁰¹²³⁴-⁹]/g, c => SUPD.indexOf(c)).replace(/⁺/g, '+').replace(/⁻/g, '-').replace(/[·•*]/g, '.');
  if (/^e-?$/.test(f)) return { atoms: {}, charge: -1, phase, formula: 'e-' };
  let charge = 0;
  const sc = splitCharge(f); f = sc.body; charge = sc.charge;
  f = f.replace(/[\[{]/g, '(').replace(/[\]}]/g, ')');
  const atoms = {}; const add = (a, k) => { for (const e in a) atoms[e] = (atoms[e] || 0) + a[e] * k; };
  for (const part of f.split('.')) {
    if (!part) continue; const pm = part.match(/^(\d+(?:\.\d+)?)?(.*)$/); const mult = pm[1] ? +pm[1] : 1; let i = 0; const s = pm[2];
    const group = () => { const a = {}; while (i < s.length && s[i] !== ')') { if (s[i] === '(') { i++; const inner = group(); i++; const n = num(); for (const e in inner) a[e] = (a[e] || 0) + inner[e] * n; } else { const em = s.slice(i).match(/^[A-Z][a-z]?/); if (!em || !ELEMENTS[em[0]]) throw new Error(`Unknown element in "${input}" near "${s.slice(i, i + 3)}"`); i += em[0].length; const n = num(); a[em[0]] = (a[em[0]] || 0) + n; } } return a; };
    const num = () => { const nm = s.slice(i).match(/^\d+(\.\d+)?/); if (nm) { i += nm[0].length; return +nm[0]; } return 1; };
    add(group(), mult);
  }
  return { atoms, charge, phase, formula: f };
}
export function molarMass(formula) { const { atoms } = typeof formula === 'string' ? parseFormula(formula) : formula; let m = 0; for (const e in atoms) m += ELEMENTS[e].mass * atoms[e]; return m; }
export function massFractions(formula) { const { atoms } = parseFormula(formula), M = molarMass({ atoms }); const o = {}; for (const e in atoms) o[e] = ELEMENTS[e].mass * atoms[e] / M; return o; }
export function prettyFormula(input) {
  const s = String(input).trim(); const ph = s.match(/\((s|l|g|aq)\)$/); let f = ph ? s.slice(0, ph.index) : s; let chg = '';
  if (/^e-?$/.test(f)) return 'e⁻';
  const sc = splitCharge(f); if (sc.charge) { const n = Math.abs(sc.charge); chg = superscript((n !== 1 ? n : '') + (sc.charge > 0 ? '+' : '-')); f = sc.body; }
  f = f.replace(/(?<=[A-Za-z)\]])(\d+(\.\d+(?![A-Za-z(]))?)/g, d => subscript(d)).replace(/[.*]/g, '·');
  return f + chg + (ph ? `(${ph[1]})` : '');
}

// ---- exact rational nullspace for balancing ----
const gcd = (a, b) => { a = Math.abs(a); b = Math.abs(b); while (b) [a, b] = [b, a % b]; return a || 1; };
const frac = (n, d = 1) => { if (d < 0) { n = -n; d = -d; } const g = gcd(n, d); return [n / g, d / g]; };
const fsub = (a, b) => frac(a[0] * b[1] - b[0] * a[1], a[1] * b[1]), fmul = (a, b) => frac(a[0] * b[0], a[1] * b[1]), fdiv = (a, b) => frac(a[0] * b[1], a[1] * b[0]);
// Balance reactants -> products. Species are formula strings. Returns integer coefficients.
export function balance(reactants, products) {
  const sp = [...reactants, ...products].map(parseFormula), nR = reactants.length, els = [...new Set(sp.flatMap(s => Object.keys(s.atoms)))];
  const SC = 1000; // support fractional stoichiometries to 3 decimals
  const rows = els.map(e => sp.map((s, j) => frac(Math.round((s.atoms[e] || 0) * SC) * (j < nR ? 1 : -1), 1)));
  if (sp.some(s => s.charge)) rows.push(sp.map((s, j) => frac(s.charge * SC * (j < nR ? 1 : -1), 1)));
  const n = sp.length, piv = []; let r = 0;
  for (let c = 0; c < n && r < rows.length; c++) {
    let p = rows.findIndex((row, i) => i >= r && row[c][0] !== 0); if (p < 0) continue;[rows[r], rows[p]] = [rows[p], rows[r]];
    const pv = rows[r][c]; rows[r] = rows[r].map(v => fdiv(v, pv));
    for (let i = 0; i < rows.length; i++) if (i !== r && rows[i][c][0] !== 0) { const f = rows[i][c]; rows[i] = rows[i].map((v, j) => fsub(v, fmul(f, rows[r][j]))); }
    piv.push(c); r++;
  }
  const free = [...Array(n).keys()].filter(c => !piv.includes(c));
  if (free.length !== 1) throw new Error(free.length === 0 ? 'Equation cannot be balanced (no non-trivial solution): check species.' : 'Equation is under-determined (multiple independent reactions): split it into separate steps.');
  const x = Array(n).fill(null); x[free[0]] = [1, 1]; piv.forEach((c, i) => { x[c] = frac(-rows[i][free[0]][0], rows[i][free[0]][1]); });
  const lcm = x.reduce((l, v) => l * v[1] / gcd(l, v[1]), 1); let ints = x.map(v => v[0] * lcm / v[1]); const g = ints.reduce((a, b) => gcd(a, b), 0); ints = ints.map(v => v / g);
  if (ints.some(v => v < 0)) { if (ints.every(v => v <= 0)) ints = ints.map(v => -v); else throw new Error('Balancing gives negative coefficients: a species is on the wrong side.'); }
  if (ints.some(v => v === 0)) throw new Error('A species has zero coefficient: it does not take part in the reaction.');
  return { coefficients: ints, reactants: ints.slice(0, nR), products: ints.slice(nR) };
}
const term = (c, s) => (c === 1 ? '' : c + ' ') + prettyFormula(s);
export function equationString(reactants, products, coeffs, arrow = '→') {
  const c = coeffs || balance(reactants, products).coefficients, nR = reactants.length;
  return reactants.map((s, i) => term(c[i], s)).join(' + ') + ` ${arrow} ` + products.map((s, i) => term(c[nR + i], s)).join(' + ');
}

// Strong electrolytes that dissociate fully in water: formula -> [[ion, count], ...]
export const DISSOCIATION = {
  'Zn(NO3)2': [['Zn2+', 1], ['NO3-', 2]], 'ZnCl2': [['Zn2+', 1], ['Cl-', 2]], 'ZnSO4': [['Zn2+', 1], ['SO4^2-', 1]], 'Zn(CH3COO)2': [['Zn2+', 1], ['CH3COO-', 2]],
  'NaOH': [['Na+', 1], ['OH-', 1]], 'KOH': [['K+', 1], ['OH-', 1]], 'LiOH': [['Li+', 1], ['OH-', 1]], 'NH4OH': [['NH4+', 1], ['OH-', 1]],
  'NaCl': [['Na+', 1], ['Cl-', 1]], 'KCl': [['K+', 1], ['Cl-', 1]], 'NaNO3': [['Na+', 1], ['NO3-', 1]], 'KNO3': [['K+', 1], ['NO3-', 1]], 'Na2SO4': [['Na+', 2], ['SO4^2-', 1]], 'NH4NO3': [['NH4+', 1], ['NO3-', 1]], 'NH4Cl': [['NH4+', 1], ['Cl-', 1]],
  'HCl': [['H+', 1], ['Cl-', 1]], 'HNO3': [['H+', 1], ['NO3-', 1]], 'H2SO4': [['H+', 2], ['SO4^2-', 1]],
  'AgNO3': [['Ag+', 1], ['NO3-', 1]], 'HAuCl4': [['H+', 1], ['AuCl4-', 1]], 'CuSO4': [['Cu2+', 1], ['SO4^2-', 1]], 'Cu(NO3)2': [['Cu2+', 1], ['NO3-', 2]],
  'FeCl3': [['Fe3+', 1], ['Cl-', 3]], 'FeCl2': [['Fe2+', 1], ['Cl-', 2]], 'FeSO4': [['Fe2+', 1], ['SO4^2-', 1]], 'Fe(NO3)3': [['Fe3+', 1], ['NO3-', 3]],
  'NaBH4': [['Na+', 1], ['BH4-', 1]], 'Na3C6H5O7': [['Na+', 3], ['C6H5O7^3-', 1]], 'CdCl2': [['Cd2+', 1], ['Cl-', 2]], 'Na2S': [['Na+', 2], ['S2-', 1]], 'Na2SeO3': [['Na+', 2], ['SeO3^2-', 1]],
  'TiCl4': [['Ti4+', 1], ['Cl-', 4]], 'Ce(NO3)3': [['Ce3+', 1], ['NO3-', 3]], 'Ni(NO3)2': [['Ni2+', 1], ['NO3-', 2]], 'Co(NO3)2': [['Co2+', 1], ['NO3-', 2]], 'NiCl2': [['Ni2+', 1], ['Cl-', 2]], 'CoCl2': [['Co2+', 1], ['Cl-', 2]],
  'Al(NO3)3': [['Al3+', 1], ['NO3-', 3]], 'AlCl3': [['Al3+', 1], ['Cl-', 3]], 'MgCl2': [['Mg2+', 1], ['Cl-', 2]], 'CaCl2': [['Ca2+', 1], ['Cl-', 2]], 'Na2CO3': [['Na+', 2], ['CO3^2-', 1]], 'BaCl2': [['Ba2+', 1], ['Cl-', 2]],
  'SnCl2': [['Sn2+', 1], ['Cl-', 2]], 'SnCl4': [['Sn4+', 1], ['Cl-', 4]], 'In(NO3)3': [['In3+', 1], ['NO3-', 3]], 'Na2WO4': [['Na+', 2], ['WO4^2-', 1]], 'Na2MoO4': [['Na+', 2], ['MoO4^2-', 1]], 'K2PtCl4': [['K+', 2], ['PtCl4^2-', 1]], 'PdCl2': [['Pd2+', 1], ['Cl-', 2]],
};
// Build molecular, complete-ionic and net-ionic equations. species: [{f:'Zn(NO3)2', phase:'aq'|'s'|'l'|'g'}]
export function ionicEquations(reactants, products) {
  const rf = reactants.map(s => s.f), pf = products.map(s => s.f), bal = balance(rf, pf), nR = rf.length;
  const withPh = (s) => s.f + (s.phase ? `(${s.phase})` : '');
  const molecular = reactants.map((s, i) => term(bal.coefficients[i], withPh(s))).join(' + ') + ' → ' + products.map((s, i) => term(bal.coefficients[nR + i], withPh(s))).join(' + ');
  const expand = (list, coefs) => { const m = new Map(); list.forEach((s, i) => { const d = s.phase === 'aq' && DISSOCIATION[s.f]; if (d) for (const [ion, k] of d) m.set(ion + '(aq)', (m.get(ion + '(aq)') || 0) + k * coefs[i]); else m.set(withPh(s), (m.get(withPh(s)) || 0) + coefs[i]); }); return m; };
  const L = expand(reactants, bal.reactants), Rr = expand(products, bal.products);
  const show = m => [...m].filter(([, c]) => c > 0).map(([s, c]) => term(c, s)).join(' + ');
  const completeIonic = show(L) + ' → ' + show(Rr);
  const spectators = [];
  for (const [s, c] of L) if (Rr.has(s)) { const k = Math.min(c, Rr.get(s)); spectators.push(prettyFormula(s)); L.set(s, c - k); Rr.set(s, Rr.get(s) - k); }
  let all = [...L.values(), ...Rr.values()].filter(v => v > 0); const g = all.reduce((a, b) => gcd(a, b), 0) || 1;
  for (const m of [L, Rr]) for (const [s, c] of m) m.set(s, c / g);
  const netIonic = all.length ? show(L) + ' → ' + show(Rr) : '(no net reaction: all species are spectators)';
  return { molecular, completeIonic, netIonic, spectators, coefficients: bal.coefficients };
}
// Charge & mass audit of a written equation
export function auditEquation(reactants, products, coeffs) {
  const sp = [...reactants, ...products].map(parseFormula), nR = reactants.length, tot = {}; let q = 0;
  sp.forEach((s, j) => { const k = coeffs[j] * (j < nR ? 1 : -1); for (const e in s.atoms) tot[e] = (tot[e] || 0) + s.atoms[e] * k; q += s.charge * k; });
  const bad = Object.entries(tot).filter(([, v]) => Math.abs(v) > 1e-9).map(([e]) => e);
  return { massBalanced: bad.length === 0, chargeBalanced: Math.abs(q) < 1e-9, unbalancedElements: bad, netCharge: q };
}
// Stoichiometric yield helper
export function theoreticalYield(reactants, products, molesAvailable, productIndex = 0) {
  const b = balance(reactants, products); let xi = Infinity, limiting = 0;
  b.reactants.forEach((c, i) => { const e = molesAvailable[i] / c; if (e < xi) { xi = e; limiting = i; } });
  const n = xi * b.products[productIndex]; return { extent: xi, limiting: reactants[limiting], moles: n, mass_g: n * molarMass(products[productIndex]), coefficients: b.coefficients };
}
