// Self-test for js/geometry/atomic.js: crystal data, builders, analysis, diffraction, writers and parsers.
import * as A from '../js/geometry/atomic.js';

export default async function () {
  const log = []; let pass = 0, fail = 0;
  const ok = (cond, name, detail = '') => { if (cond) pass++; else { fail++; log.push(`✗ ${name}${detail ? ' — ' + detail : ''}`); } };
  const near = (v, ref, tol, name) => ok(Number.isFinite(v) && Math.abs(v - ref) <= tol, name, `got ${v}, expected ${ref} ± ${tol}`);
  const throws = (fn, re, name) => { try { fn(); ok(false, name, 'did not throw'); } catch (e) { ok(e instanceof Error && (!re || re.test(e.message)), name, `message: ${e.message}`); } };
  const time = fn => { let b = Infinity, r; for (let i = 0; i < 5; i++) { const t = performance.now(); r = fn(); b = Math.min(b, performance.now() - t); } return [b, r]; };

  // ── exports present
  const EXPORTS = ['CRYSTALS', 'latticeVectors', 'buildCell', 'supercell', 'nanoparticle', 'nanosheet', 'nanotube', 'coreShell', 'addDefects', 'slab', 'analyse', 'dSpacings', 'powderXRD', 'rdf',
    'toCIF', 'toXYZ', 'toExtXYZ', 'toPOSCAR', 'toXSF', 'toPDB', 'toLAMMPSData', 'toMOL',
    'parseCIF', 'parseXYZ', 'parsePOSCAR', 'parseXSF', 'parsePDB', 'parseMOL', 'parseSDF', 'parseMOL2', 'parseLAMMPSData', 'parseGRO', 'parseCML', 'parseGaussianInput', 'parseQEInput', 'parseCASTEPCell', 'parseLAMMPSDump'];
  for (const e of EXPORTS) ok(A[e] !== undefined, `export ${e}`);
  const KEYS = ['ZnO_wurtzite', 'ZnO_zincblende', 'TiO2_rutile', 'TiO2_anatase', 'Si', 'Ge', 'C_diamond', 'graphite', 'graphene', 'hBN', 'MoS2_2H', 'GaAs', 'GaN_wurtzite', 'CdSe_wurtzite', 'CdS_zincblende', 'ZnS_zincblende', 'Au', 'Ag', 'Cu', 'Pt', 'Pd', 'Ni', 'Al', 'Fe_bcc', 'W', 'NaCl', 'MgO', 'CeO2', 'SrTiO3', 'BaTiO3_cubic', 'CsCl', 'Fe3O4', 'SiO2_cristobalite', 'Al2O3_corundum'];
  for (const k of KEYS) ok(!!A.CRYSTALS[k], `crystal ${k}`);

  // ── every library crystal: stoichiometry, sane contacts, serialisable
  const an = {};
  for (const k of Object.keys(A.CRYSTALS)) {
    const d = A.CRYSTALS[k], s = A.buildCell(k), a = an[k] = A.analyse(s);
    ok(a.formula === d.formula, `${k} cell stoichiometry`, `${a.formula} vs ${d.formula}`);
    ok(a.nearestNeighbour_A > 1.3 && a.nearestNeighbour_A < 3.7, `${k} nearest-neighbour distance plausible`, String(a.nearestNeighbour_A));
    ok(a.surfaceFraction === 0, `${k} bulk has no under-coordinated atoms`, String(a.surfaceFraction));
    ok(d.basis.every(b => b.f.every(v => v >= 0 && v < 1)), `${k} basis inside the cell`);
    ok(JSON.parse(JSON.stringify(s)).atoms.length === s.atoms.length, `${k} JSON-serialisable`);
  }
  ok(A.CRYSTALS.Fe3O4.basis.length === 56, 'Fe3O4 has 56 atoms'); ok(A.CRYSTALS.Al2O3_corundum.basis.length === 30, 'corundum has 30 atoms');
  ok(A.CRYSTALS.TiO2_anatase.basis.length === 12 && A.CRYSTALS.SiO2_cristobalite.basis.length === 24, 'anatase 12 / cristobalite 24 atoms');

  // ── densities against handbook values (1.5 %)
  for (const [k, rho] of [['ZnO_wurtzite', 5606], ['Si', 2329], ['Au', 19300], ['TiO2_rutile', 4250], ['NaCl', 2165], ['Cu', 8960], ['Al', 2700], ['Fe_bcc', 7874], ['MgO', 3580], ['Al2O3_corundum', 3987], ['Fe3O4', 5200], ['TiO2_anatase', 3890], ['C_diamond', 3515], ['graphite', 2260], ['CeO2', 7215], ['GaAs', 5320], ['W', 19250], ['Ag', 10490], ['SrTiO3', 5120], ['MoS2_2H', 5060]])
    ok(Math.abs(an[k].density_kg_m3 / rho - 1) < 0.015, `${k} density within 1.5 % of ${rho}`, an[k].density_kg_m3.toFixed(0));

  // ── bond lengths and coordination
  const bond = (k, pair) => an[k].bonds.find(b => b.pair === pair);
  const zb = bond('ZnO_wurtzite', 'Zn-O');
  ok(zb && zb.min_A > 1.97 && zb.max_A < 1.995, 'ZnO Zn–O bonds 1.97–1.99 Å', JSON.stringify(zb)); ok(an.ZnO_wurtzite.coordination.Zn === 4 && an.ZnO_wurtzite.coordination.O === 4, 'ZnO tetrahedral');
  near(bond('TiO2_rutile', 'Ti-O').min_A, 1.949, 0.006, 'rutile equatorial Ti–O'); near(bond('TiO2_rutile', 'Ti-O').max_A, 1.980, 0.006, 'rutile apical Ti–O');
  near(bond('TiO2_anatase', 'Ti-O').min_A, 1.934, 0.006, 'anatase equatorial Ti–O'); near(bond('TiO2_anatase', 'Ti-O').max_A, 1.980, 0.006, 'anatase apical Ti–O');
  ok(an.TiO2_rutile.coordination.Ti === 6 && an.TiO2_rutile.coordination.O === 3 && an.TiO2_anatase.coordination.Ti === 6, 'TiO2 6:3 coordination');
  near(bond('Al2O3_corundum', 'Al-O').min_A, 1.855, 0.008, 'corundum short Al–O'); near(bond('Al2O3_corundum', 'Al-O').max_A, 1.972, 0.008, 'corundum long Al–O');
  ok(an.Al2O3_corundum.coordination.Al === 6 && an.Al2O3_corundum.coordination.O === 4, 'corundum 6:4 coordination');
  near(bond('Fe3O4', 'Fe-O').min_A, 1.889, 0.01, 'magnetite tetrahedral Fe–O'); near(bond('Fe3O4', 'Fe-O').max_A, 2.059, 0.01, 'magnetite octahedral Fe–O');
  near(an.Fe3O4.coordination.Fe, (8 * 4 + 16 * 6) / 24, 1e-9, 'magnetite mean Fe coordination'); ok(an.Fe3O4.coordination.O === 4, 'magnetite O coordination');
  near(bond('Si', 'Si-Si').length_A, 2.352, 0.002, 'Si–Si bond'); near(bond('Au', 'Au-Au').length_A, 2.884, 0.002, 'Au–Au distance'); near(bond('NaCl', 'Na-Cl').length_A, 2.820, 0.002, 'Na–Cl distance');
  near(bond('MoS2_2H', 'Mo-S').length_A, 2.41, 0.02, 'Mo–S bond'); near(bond('graphite', 'C-C').length_A, 1.422, 0.003, 'graphite C–C');
  ok(an.Au.coordination.Au === 12 && an.Fe_bcc.coordination.Fe === 8 && an.Si.coordination.Si === 4 && an.NaCl.coordination.Na === 6 && an.CsCl.coordination.Cs === 8 && an.CeO2.coordination.Ce === 8 && an.graphene.coordination.C === 3 && an.MoS2_2H.coordination.Mo === 6 && an.SrTiO3.coordination.Ti === 6 && an.SrTiO3.coordination.Sr === 12, 'coordination numbers of reference structures');

  // ── lattice vectors, supercell, definitions
  const lv = A.latticeVectors({ a: 3, b: 4, c: 5, alpha: 80, beta: 95, gamma: 110 }), cp = A.cellParams(lv);
  ok([cp.a - 3, cp.b - 4, cp.c - 5, cp.alpha - 80, cp.beta - 95, cp.gamma - 110].every(v => Math.abs(v) < 1e-9), 'latticeVectors ↔ cellParams triclinic round trip');
  throws(() => A.latticeVectors({ a: 3, b: 3, c: 3, alpha: 10, beta: 100, gamma: 140 }), /valid cell/, 'impossible cell angles rejected');
  throws(() => A.buildCell('unobtainium'), /Unknown crystal/, 'unknown crystal key rejected');
  const zno = A.buildCell('ZnO_wurtzite'), sc = A.supercell(zno, [3, 2, 2]), asc = A.analyse(sc);
  ok(sc.atoms.length === 48, 'supercell atom count'); near(asc.density_kg_m3, an.ZnO_wurtzite.density_kg_m3, 1e-6, 'supercell keeps density'); ok(asc.coordination.Zn === 4, 'supercell keeps coordination');
  const custom = A.buildCell({ name: 'CsCl-type test', a: 4, basis: [{ el: 'Cs', f: [0, 0, 0] }, { el: 'Cl', f: [0.5, 0.5, 0.5] }] });
  ok(custom.atoms.length === 2 && custom.formula === 'CsCl', 'custom definition builds');

  // ── nanoparticles
  const [tNp, np] = time(() => A.nanoparticle('ZnO_wurtzite', { shape: 'sphere', diameter_nm: 5 }));
  ok(tNp < 500, '5 nm nanoparticle carved in < 500 ms', `${tNp.toFixed(0)} ms`);
  near(np.atoms.length, Math.PI / 6 * 50 ** 3 * 4 / 47.61, 150, '5 nm ZnO sphere atom count matches bulk number density');
  ok(np.cell === null && np.truncated === false && np.atoms.every(a => Math.hypot(a.x, a.y, a.z) <= 25.001), 'sphere atoms inside radius');
  const anp = A.analyse(np); ok(anp.surfaceFraction > 0.08 && anp.surfaceFraction < 0.4, 'nanoparticle surface fraction plausible', String(anp.surfaceFraction)); near(anp.radiusOfGyration_A, Math.sqrt(3 / 5) * 25, 0.6, 'sphere radius of gyration = √(3/5) R'); ok(anp.density_kg_m3 === null, 'finite cluster has no cell density');
  for (const shape of ['cube', 'rod', 'hexrod', 'octahedron', 'wulff-like']) for (const k of ['Au', 'ZnO_wurtzite']) { const p = A.nanoparticle(k, { shape, diameter_nm: 2.5, length_nm: 5 }); ok(p.atoms.length > 50 && A.analyse(p).nearestNeighbour_A > 1.9, `${k} ${shape} particle`, String(p.atoms.length)); }
  const hex = A.nanoparticle('ZnO_wurtzite', { shape: 'hexrod', diameter_nm: 3, length_nm: 8 });
  ok(hex.atoms.every(a => Math.abs(a.z) <= 40.001 && Math.abs(a.y) <= 15.001) && Math.max(...hex.atoms.map(a => a.z)) - Math.min(...hex.atoms.map(a => a.z)) > 74, 'hexagonal nanorod along c with requested length');
  const [tBig, big] = time(() => A.nanoparticle('Au', { diameter_nm: 60 }));
  ok(big.truncated === true && big.atoms.length <= A.MAX_ATOMS && big.diameter_nm < 60 && tBig < 1500, 'oversized particle is capped and flagged truncated', `${big.atoms.length} atoms, ${tBig.toFixed(0)} ms`);
  throws(() => A.nanoparticle('Au', { shape: 'banana' }), /unknown shape/, 'unknown shape rejected');

  // ── nanotubes
  const nt = A.nanotube({ n: 10, m: 10, length_nm: 3 }), ant = A.analyse(nt);
  near(nt.diameter_nm, 1.356, 0.001, '(10,10) nanotube diameter'); near(nt.chiralAngle_deg, 30, 1e-6, 'armchair chiral angle'); ok(nt.metallic === true, 'armchair tube metallic');
  ok(ant.coordination.C === 3 && nt.atoms.length % 40 === 0, 'periodic armchair tube fully three-coordinated'); near(ant.bonds[0].length_A, 1.42, 0.005, 'nanotube C–C bond');
  const zz = A.nanotube({ n: 10, m: 0, length_nm: 2 }); ok(zz.metallic === false && Math.abs(zz.chiralAngle_deg) < 1e-9, '(10,0) zigzag semiconducting'); near(zz.diameter_nm, 0.783, 0.001, '(10,0) diameter');
  const ch = A.nanotube({ n: 8, m: 4, length_nm: 4 }); ok(A.analyse(ch).coordination.C === 3 && ch.atomsPerUnitCell === 112, '(8,4) chiral tube periodic and three-coordinated');
  const bnt = A.nanotube({ n: 6, m: 6, length_nm: 2, element: 'BN', bond_A: 1.45 }); ok(A.analyse(bnt).formula === 'BN' && bnt.metallic === false, 'BN nanotube');
  throws(() => A.nanotube({ n: 0, m: 0 }), /chiral indices/, 'invalid chiral indices rejected');

  // ── slabs and sheets
  const s111 = A.slab('Au', { hkl: [1, 1, 1], layers: 3, vacuum_A: 12 }); near(s111.d_hkl_A, 4.0782 / Math.sqrt(3), 1e-3, 'Au(111) interplanar spacing'); ok(s111.atoms.length === 12, 'Au(111) 3-layer slab atom count'); near(s111.cell[2][2], 2 * 4.0782 / Math.sqrt(3) + 12, 1e-3, 'slab cell = thickness + vacuum');
  const s001 = A.slab('Au', { hkl: [0, 0, 1], layers: 4 }); near(s001.thickness_A, 1.5 * 4.0782, 1e-3, 'Au(001) 4 atomic layers'); near(A.analyse(s001).nearestNeighbour_A, 2.884, 0.002, 'slab keeps bond length');
  const s110 = A.slab('Au', { hkl: [1, 1, 0], layers: 2 }); near(s110.thickness_A, 4.0782 / (2 * Math.sqrt(2)), 1e-3, 'Au(110) layer spacing');
  const s100 = A.slab('Si', { hkl: [1, 0, 0], layers: 4 }); near(s100.thickness_A, 3 * 5.431 / 4, 1e-3, 'Si(100) atomic layer spacing a/4');
  const sm = A.slab('ZnO_wurtzite', { hkl: [1, 0, -1, 0], layers: 4 }); near(sm.d_hkl_A, 3.2495 * Math.sqrt(3) / 2, 1e-3, 'ZnO (10-10) spacing'); ok(A.analyse(sm).formula === 'ZnO', 'ZnO m-plane slab stoichiometric');
  const sc0 = A.slab('ZnO_wurtzite', { hkl: [0, 0, 0, 1], layers: 8 }); near(A.analyse(sc0).bonds[0].min_A, 1.974, 0.003, 'ZnO (0001) slab bonds');
  ok(A.analyse(s111).surfaceFraction > 0.5, 'slab surface atoms detected');
  const mos = A.nanosheet('MoS2_2H', { nx: 3, ny: 3, layers: 1 }); ok(mos.atoms.length === 27 && A.analyse(mos).formula === 'MoS2' && A.analyse(mos).coordination.Mo === 6, 'MoS2 monolayer is one S–Mo–S trilayer');
  const gr = A.nanosheet('graphite', { nx: 4, ny: 4, layers: 2 }); ok(gr.atoms.length === 64, 'bilayer graphene from graphite'); near(Math.max(...gr.atoms.map(a => a.z)) - Math.min(...gr.atoms.map(a => a.z)), 3.3555, 1e-3, 'bilayer spacing');
  const g1 = A.nanosheet('graphene', { nx: 5, ny: 5 }); ok(g1.atoms.length === 50 && A.analyse(g1).coordination.C === 3, 'graphene sheet');
  throws(() => A.slab('Au', { hkl: [0, 0, 0] }), /hkl/, 'zero Miller indices rejected');

  // ── core–shell and defects
  const cs = A.coreShell('CdSe_wurtzite', 'ZnS_zincblende', { core_nm: 3, shell_nm: 0.6 });
  ok(cs.nCore > 300 && cs.nShell > 300 && cs.atoms.filter(a => a.region === 'core').every(a => Math.hypot(a.x, a.y, a.z) <= 15.001) && cs.atoms.every(a => Math.hypot(a.x, a.y, a.z) <= 21.001), 'core–shell regions'); ok(A.analyse(cs).nearestNeighbour_A > 1.75, 'core–shell interface has no overlapping atoms', String(A.analyse(cs).nearestNeighbour_A));
  const big0 = A.supercell(zno, [5, 5, 4]), df = A.addDefects(big0, { vacancyFraction: 0.05, substitution: { from: 'Zn', to: 'Al', fraction: 0.1 }, seed: 7 }), df2 = A.addDefects(big0, { vacancyFraction: 0.05, substitution: { from: 'Zn', to: 'Al', fraction: 0.1 }, seed: 7 }), df3 = A.addDefects(big0, { vacancyFraction: 0.05, seed: 8 });
  ok(df.atoms.length === 380 && df.defects.vacancies === 20 && df.atoms.filter(a => a.el === 'Al').length === df.defects.substitutions && df.defects.substitutions > 0, 'defect counts'); ok(JSON.stringify(df.atoms) === JSON.stringify(df2.atoms), 'defects reproducible for a given seed'); ok(JSON.stringify(df.atoms) !== JSON.stringify(df3.atoms) && big0.atoms.length === 400, 'different seed differs; input untouched');

  // ── d-spacings and powder XRD
  const ds = A.dSpacings('ZnO_wurtzite', { dmin_A: 1.3 }); near(ds[0].d_A, 2.8141, 0.001, 'ZnO d(100)'); ok(ds[0].hkl.join() === '1,0,0' && ds[0].multiplicity === 6 && ds[1].hkl.join() === '0,0,2' && ds[2].multiplicity === 12, 'ZnO reflection order and multiplicities');
  const dsi = A.dSpacings('Si', { dmin_A: 1.0 }); ok(dsi.map(d => d.hkl.join('')).join(' ') === '111 220 311 400 331 422 511', 'diamond-lattice extinctions (no 200, no 222)', dsi.map(d => d.hkl.join('')).join(' ')); ok(dsi[0].multiplicity === 8 && dsi[1].multiplicity === 12, 'cubic multiplicities');
  ok(A.dSpacings('Si', { dmin_A: 2, includeForbidden: true }).some(d => d.hkl.join('') === '100' && d.allowed === false), 'forbidden reflections available on request');
  const [tX, xz] = time(() => A.powderXRD('ZnO_wurtzite', { twoThetaMin: 10, twoThetaMax: 90 }));
  ok(tX < 300, 'ZnO powder pattern in < 300 ms', `${tX.toFixed(0)} ms`); ok(xz.twoTheta.length === 4001 && xz.intensity.length === 4001 && xz.intensity.every(Number.isFinite) && Math.abs(Math.max(...xz.intensity) - 100) < 1e-6, 'pattern grid and normalisation');
  const top3 = [...xz.peaks].sort((a, b) => b.I_rel - a.I_rel).slice(0, 3).map(p => p.twoTheta).sort((a, b) => a - b);
  near(top3[0], 31.77, 0.1, 'ZnO (100) at 31.8°'); near(top3[1], 34.42, 0.1, 'ZnO (002) at 34.4°'); near(top3[2], 36.25, 0.1, 'ZnO (101) at 36.3°');
  const pz = h => xz.peaks.find(p => p.hkl.join('') === h); ok(pz('101').I_rel === 100 && pz('100').I_rel > 45 && pz('100').I_rel < 70 && pz('002').I_rel > 30 && pz('002').I_rel < 55, 'ZnO relative intensities (reference ≈ 57 : 44 : 100)', `${pz('100').I_rel.toFixed(0)} : ${pz('002').I_rel.toFixed(0)} : 100`);
  near(xz.twoTheta[xz.intensity.indexOf(Math.max(...xz.intensity))], 36.25, 0.06, 'profile maximum at (101)');
  const xs = A.powderXRD('Si'); ok(xs.peaks[0].hkl.join('') === '111' && xs.peaks[0].I_rel === 100, 'Si strongest peak is (111)'); near(xs.peaks[0].twoTheta, 28.44, 0.02, 'Si (111) at 28.44°'); ok(xs.peaks.map(p => p.hkl.join('')).join(' ') === '111 220 311 400 331 422', 'Si peak list 10–90°');
  const xa = A.powderXRD('Au'); near(xa.peaks[0].twoTheta, 38.19, 0.03, 'Au (111) at 38.2°'); ok(xa.peaks.length === 5 && xa.peaks[1].I_rel > 35 && xa.peaks[1].I_rel < 60, 'Au fcc peaks and (200) intensity');
  const xn = A.powderXRD('NaCl'); ok(xn.peaks.find(p => p.I_rel === 100).hkl.join('') === '200', 'NaCl strongest peak is (200)');
  const xan = A.powderXRD('TiO2_anatase'); near(xan.peaks[0].twoTheta, 25.31, 0.05, 'anatase (101) at 25.3°'); ok(xan.peaks[0].I_rel === 100, 'anatase (101) strongest');
  const xr = A.powderXRD('TiO2_rutile'); near(xr.peaks.find(p => p.I_rel === 100).twoTheta, 27.44, 0.05, 'rutile (110) at 27.4° strongest');
  const xf = A.powderXRD('Fe3O4'); near(xf.peaks.find(p => p.I_rel === 100).twoTheta, 35.43, 0.05, 'magnetite (311) at 35.4° strongest');
  const fwhm = (x, c) => { const i0 = x.twoTheta.findIndex(t => t >= c - 3), i1 = x.twoTheta.findIndex(t => t >= c + 3); let m = 0, im = i0; for (let i = i0; i < i1; i++) if (x.intensity[i] > m) { m = x.intensity[i]; im = i; } let l = im, r = im; while (x.intensity[l] > m / 2) l--; while (x.intensity[r] > m / 2) r++; return x.twoTheta[r] - x.twoTheta[l]; };
  const xb = A.powderXRD('Si', { crystallite_nm: 10, instrument_deg: 0 }), scherrer = 0.9 * 1.5406 / (100 * Math.cos(14.22 * Math.PI / 180)) * 180 / Math.PI;
  near(fwhm(xb, 28.44), scherrer, 0.06, 'Scherrer width for 10 nm crystallites'); ok(fwhm(A.powderXRD('Si', { crystallite_nm: 10, strain: 0.01 }), 28.44) > fwhm(xb, 28.44) + 0.2, 'microstrain broadens the peak');
  near(A.powderXRD('Si', { wavelength_A: 0.7093 }).peaks[0].twoTheta, 2 * Math.asin(0.7093 / (2 * 5.431 / Math.sqrt(3))) * 180 / Math.PI, 0.01, 'Mo Kα shifts Si (111)');
  for (const el of ['H', 'C', 'O', 'Si', 'Ti', 'Fe', 'Zn', 'Mo', 'Ag', 'Ce', 'W', 'Pt', 'Au', 'Pb']) { const Z = { H: 1, C: 6, O: 8, Si: 14, Ti: 22, Fe: 26, Zn: 30, Mo: 42, Ag: 47, Ce: 58, W: 74, Pt: 78, Au: 79, Pb: 82 }[el]; near(A.scatteringFactor(el, 0), Z, 0.06, `f(${el}, s=0) = Z`); ok(A.scatteringFactor(el, 0.5) < A.scatteringFactor(el, 0.2) && A.scatteringFactor(el, 0.2) < Z, `f(${el}) decreases with s`); }
  near(A.scatteringFactor('Hf', 0), 72, 1e-6, 'Thomas–Fermi fallback f(0) = Z'); ok(Math.abs(A.scatteringFactor('Zr', 0.3) / (A.scatteringFactor('Hf', 0.3) * 40 / 72) - 1) < 0.35, 'tabulated and Thomas–Fermi factors agree in order of magnitude');
  throws(() => A.powderXRD(np), /periodic/, 'XRD of a finite cluster rejected with a clear message');

  // ── RDF
  const g = A.rdf(A.buildCell('Si'), { rmax_A: 6, dr_A: 0.05 }), firstPk = g.r[g.g.findIndex(v => v > 1)];
  near(firstPk, 2.352, 0.05, 'Si RDF first shell'); near(g.cn[g.r.findIndex(r => r > 3)], 4, 1e-9, 'Si RDF first-shell coordination = 4'); near(g.cn[g.r.findIndex(r => r > 4.2)], 16, 1e-9, 'Si RDF second shell (12)');
  const gn = A.rdf(np, { rmax_A: 8 }); ok(gn.g.every(Number.isFinite) && /cluster/.test(gn.normalisation), 'cluster RDF finite');

  // ── writer → parser round trips
  const same = (s0, s1, name, { tol = 2e-3, cell = true } = {}) => {
    if (s1.atoms.length !== s0.atoms.length) return ok(false, name, `atom count ${s1.atoms.length} vs ${s0.atoms.length}`);
    if (cell) { if (!s1.cell) return ok(false, name, 'cell lost'); const p0 = A.cellParams(s0.cell), p1 = A.cellParams(s1.cell); for (const k of ['a', 'b', 'c', 'alpha', 'beta', 'gamma']) if (Math.abs(p0[k] - p1[k]) > 2e-3) return ok(false, name, `cell ${k}: ${p1[k]} vs ${p0[k]}`); }
    const fr = s => { if (!cell) return s.atoms.map(a => [a.el, a.x, a.y, a.z]); const M = s.cell, d = M[0][0] * (M[1][1] * M[2][2] - M[1][2] * M[2][1]) - M[0][1] * (M[1][0] * M[2][2] - M[1][2] * M[2][0]) + M[0][2] * (M[1][0] * M[2][1] - M[1][1] * M[2][0]);
      const cr = (u, v) => [u[1] * v[2] - u[2] * v[1], u[2] * v[0] - u[0] * v[2], u[0] * v[1] - u[1] * v[0]], R = [cr(M[1], M[2]), cr(M[2], M[0]), cr(M[0], M[1])].map(v => v.map(c => c / d));
      return s.atoms.map(a => [a.el, ...R.map(v => v[0] * a.x + v[1] * a.y + v[2] * a.z)]); };
    const f0 = fr(s0), f1 = fr(s1), used = new Uint8Array(f1.length), t = cell ? tol / 3 : tol;
    for (const p of f0) { let hit = -1; for (let j = 0; j < f1.length && hit < 0; j++) { if (used[j] || f1[j][0] !== p[0]) continue; let good = true; for (let k = 1; k < 4 && good; k++) { let dd = Math.abs(f1[j][k] - p[k]); if (cell) { dd %= 1; dd = Math.min(dd, 1 - dd); } if (dd > t) good = false; } if (good) hit = j; } if (hit < 0) return ok(false, name, `no match for ${p[0]} at ${p.slice(1).map(v => v.toFixed(4))}`); used[hit] = 1; }
    ok(true, name);
  };
  const mixed = A.addDefects(A.supercell(zno, [2, 2, 1]), { substitution: { from: 'Zn', to: 'Mg', fraction: 0.25 }, seed: 2 });       // 16 atoms, hexagonal cell, interleaved species
  const tri = A.buildCell({ name: 'triclinic test', a: 5.1, b: 6.2, c: 7.3, alpha: 82, beta: 97, gamma: 105, basis: [{ el: 'Cu', f: [0.1, 0.2, 0.3] }, { el: 'O', f: [0.6, 0.15, 0.8] }, { el: 'S', f: [0.35, 0.7, 0.45] }, { el: 'O', f: [0.9, 0.9, 0.05] }] });
  const cluster = A.nanoparticle('Au', { diameter_nm: 1.2 });
  const water = { kind: 'atoms', cell: null, title: 'water', atoms: [{ el: 'O', x: 0, y: 0, z: 0.1173 }, { el: 'H', x: 0, y: 0.7572, z: -0.4692 }, { el: 'H', x: 0, y: -0.7572, z: -0.4692 }] };
  for (const [nm, s] of [['ZnO:Mg supercell', mixed], ['triclinic cell', tri]]) {
    same(s, A.parseCIF(A.toCIF(s)), `CIF round trip (${nm})`); same(s, A.parsePOSCAR(A.toPOSCAR(s)), `POSCAR round trip (${nm})`); same(s, A.parseXSF(A.toXSF(s)), `XSF round trip (${nm})`);
    same(s, A.parseXYZ(A.toExtXYZ(s)), `extended XYZ round trip (${nm})`); same(s, A.parsePDB(A.toPDB(s)), `PDB round trip (${nm})`, { tol: 4e-3 });
    for (const style of ['atomic', 'charge', 'full']) same(s, A.parseLAMMPSData(A.toLAMMPSData(s, { style })), `LAMMPS data round trip, ${style} (${nm})`);
    same(s, A.parseXYZ(A.toXYZ(s)), `XYZ round trip (${nm})`, { cell: false });
  }
  for (const [nm, s] of [['Au cluster', cluster], ['water', water]]) {
    same(s, A.parseXYZ(A.toXYZ(s)), `XYZ round trip (${nm})`, { cell: false }); same(s, A.parseXYZ(A.toExtXYZ(s)), `extended XYZ round trip (${nm})`, { cell: false });
    same(s, A.parseXSF(A.toXSF(s)), `XSF ATOMS round trip (${nm})`, { cell: false }); same(s, A.parsePDB(A.toPDB(s)), `PDB round trip (${nm})`, { cell: false });
    same(s, A.parseMOL(A.toMOL(s)), `MOL round trip (${nm})`, { cell: false }); same(s, A.parseSDF(A.toMOL(s) + '$$$$\n'), `SDF round trip (${nm})`, { cell: false });
    ok(A.parseCIF(A.toCIF(s)).atoms.length === s.atoms.length && A.parsePOSCAR(A.toPOSCAR(s)).atoms.length === s.atoms.length && A.parseLAMMPSData(A.toLAMMPSData(s)).atoms.length === s.atoms.length, `boxed periodic writers accept a finite structure (${nm})`);
  }
  const wm = A.parseMOL(A.toMOL(water)); ok(wm.bonds.length === 2 && wm.bonds.every(b => b[2] === 1), 'MOL carries perceived O–H bonds'); ok(A.analyse(water).formula === 'H2O' && A.analyse(water).coordination.O === 2, 'water analysis');
  ok(A.parseXYZ(A.toExtXYZ(mixed)).format === 'extended XYZ' && A.parseLAMMPSData(A.toLAMMPSData(mixed)).atomStyle === 'atomic', 'format metadata reported');
  throws(() => A.toMOL(A.supercell(zno, [7, 7, 6])), /999/, 'MOL writer refuses > 999 atoms'); throws(() => A.toXYZ({ kind: 'atoms', atoms: [] }), /no atoms/, 'writers reject empty structures');

  // ── CIF: real-world features
  const cif = ['#------------------------------------------------------------', 'data_ZnO_zincite', "_chemical_name_mineral 'Zincite'", '_chemical_formula_sum "Zn O"', '_publ_section_title', ';', ' A multi-line text field with a loop_ keyword and data_ inside it', ';',
    '_cell_length_a\t3.2495(5)', '_cell_length_b 3.2495(5)', '_cell_length_c   5.2069(8)  # trailing comment', '_cell_angle_alpha 90', '_cell_angle_beta 90.00', '_cell_angle_gamma 120.', "_symmetry_space_group_name_H-M 'P 63 m c'", '_symmetry_Int_Tables_number 186',
    'loop_', '_symmetry_equiv_pos_site_id', '_symmetry_equiv_pos_as_xyz', "1 'x, y, z'", "2 '-y, x-y, z'", "3 '-x+y, -x, z'", "4 '-x, -y, z+1/2'", "5 'y, -x+y, 1/2+z'", '6 x-y,x,z+1/2', "7 '-y, -x, z'", "8 '-x+y, y, z'", "9 'x, x-y, z'", "10 'y, x, z+0.5'", "11 'x-y, -y, z+1/2'", "12 '-X, -X+Y, Z+1/2'",
    'loop_', '_atom_site_label', '_atom_site_type_symbol', '_atom_site_symmetry_multiplicity', '_atom_site_Wyckoff_symbol', '_atom_site_fract_x', '_atom_site_fract_y', '_atom_site_fract_z', '_atom_site_occupancy', '_atom_site_U_iso_or_equiv',
    'Zn1 Zn2+ 2 b 0.3333 0.6667 0.0000 1 0.0063(2)', 'O1 O2- 2 b 0.33333(4) 0.66667(4) 0.3819(3) 1.0 0.0068(4)', '', 'loop_', '_atom_site_aniso_label', '_atom_site_aniso_U_11', 'Zn1 0.006', 'O1 0.007', ''].join('\r\n');
  const pc = A.parseCIF(cif), apc = A.analyse(pc);
  ok(pc.atoms.length === 4 && apc.formula === 'ZnO' && pc.nSymops === 12, 'CIF with symmetry operators expands 2 sites × 12 operators to 4 atoms', `${pc.atoms.length} atoms`); near(apc.density_kg_m3, an.ZnO_wurtzite.density_kg_m3, 2, 'CIF-derived density'); near(apc.bonds[0].min_A, 1.974, 0.003, 'CIF-derived Zn–O bond'); ok(pc.title === 'Zincite' && /P 63 m c/.test(pc.spacegroup), 'CIF title and space group');
  const cifSi = "data_Si\n_cell_length_a 5.431\n_cell_length_b 5.431\n_cell_length_c 5.431\n_cell_angle_alpha 90\n_cell_angle_beta 90\n_cell_angle_gamma 90\nloop_\n_space_group_symop_operation_xyz\n'x,y,z'\n'x,y+1/2,z+1/2'\n'x+1/2,y,z+1/2'\n'x+1/2,y+1/2,z'\n'-x+1/4,-y+1/4,-z+1/4'\n'-x+1/4,-y+3/4,-z+3/4'\n'-x+3/4,-y+1/4,-z+3/4'\n'-x+3/4,-y+3/4,-z+1/4'\nloop_\n_atom_site_label\n_atom_site_fract_x\n_atom_site_fract_y\n_atom_site_fract_z\nSi1 0 0 0\n";
  const psi = A.parseCIF(cifSi); ok(psi.atoms.length === 8 && A.analyse(psi).coordination.Si === 4, 'CIF without type symbols: element from label, special positions de-duplicated');
  const mm = "data_1ABC\n_cell.length_a 20\n_cell.length_b 20\n_cell.length_c 20\n_cell.angle_alpha 90\n_cell.angle_beta 90\n_cell.angle_gamma 90\nloop_\n_atom_site.group_PDB\n_atom_site.id\n_atom_site.type_symbol\n_atom_site.label_atom_id\n_atom_site.Cartn_x\n_atom_site.Cartn_y\n_atom_site.Cartn_z\nATOM 1 N N 1.0 2.0 3.0\nATOM 2 C CA 2.4 2.0 3.0\nHETATM 3 ZN ZN 5.0 5.0 5.0\n";
  const pmm = A.parseCIF(mm); ok(pmm.atoms.length === 3 && pmm.atoms.map(a => a.el).join() === 'N,C,Zn' && pmm.atoms[1].x === 2.4, 'mmCIF-style Cartesian atom sites');
  ok(A.parseCIF("data_x\n_cell_length_a 4\n_cell_length_b 4\n_cell_length_c 4\n_symmetry_space_group_name_H-M 'F m -3 m'\nloop_\n_atom_site_label\n_atom_site_fract_x\n_atom_site_fract_y\n_atom_site_fract_z\nCu 0 0 0\n").warnings.some(w => /no symmetry operators/.test(w)), 'CIF without operators warns that the cell may be incomplete');
  throws(() => A.parseCIF('data_x\n_cell_length_a 4\n'), /_atom_site/, 'CIF without atoms rejected');
  throws(() => A.parseCIF('data_x\nloop_\n_atom_site_label\n_atom_site_fract_x\n_atom_site_fract_y\n_atom_site_fract_z\nCu 0 0 0\n'), /cell/i, 'CIF fractional coordinates without a cell rejected');
  throws(() => A.parseCIF(cifSi.replace('Si1 0 0 0', 'Si1 0 0')), /multiple|incomplete/, 'CIF with an incomplete atom row rejected');
  throws(() => A.parseCIF(cifSi.replace('Si1 0 0 0', 'Si1 0 abc 0')), /non-numeric/, 'CIF with a non-numeric coordinate rejected');
  throws(() => A.parseCIF(cifSi.replace("'x,y,z'", "'x,y'")), /symmetry operator/, 'CIF with a malformed symmetry operator rejected');
  throws(() => A.parseCIF("data_x\n_a 'unterminated\n"), /unterminated/, 'CIF with an unterminated string rejected');

  // ── other parsers on hand-written real-world style input
  const xyz2 = A.parseXYZ('3\r\nwater molecule\r\nO\t0.0\t0.0\t0.1173\r\n1 0.0 0.7572 -0.4692\r\nh 0.0 -0.7572 -4.692E-1\r\n\r\n');
  ok(xyz2.atoms.map(a => a.el).join('') === 'OHH' && xyz2.title === 'water molecule' && Math.abs(xyz2.atoms[2].z + 0.4692) < 1e-12, 'XYZ with CRLF, tabs, atomic numbers and exponents');
  const xyzT = A.parseXYZ('1\nframe 0\nAr 0 0 0\n1\nframe 1\nAr 0 0 1\n'); ok(xyzT.nFrames === 2 && xyzT.frames[1].atoms[0].z === 1, 'multi-frame XYZ');
  const exyz = A.parseXYZ('2\nLattice="4.0 0.0 0.0 0.0 4.0 0.0 0.0 0.0 4.0" Properties=id:I:1:species:S:1:pos:R:3:forces:R:3 pbc="T T T"\n1 Cs 0 0 0 0.1 0.1 0.1\n2 Cl 2 2 2 0 0 0\n'); ok(exyz.cell[1][1] === 4 && exyz.atoms[1].el === 'Cl' && exyz.atoms[1].x === 2, 'extended XYZ with extra property columns');
  throws(() => A.parseXYZ('5\ncomment\nH 0 0 0\n'), /declares 5 atoms/, 'truncated XYZ rejected'); throws(() => A.parseXYZ('2\nc\nH 0 0 0\nH 0 zero 1\n'), /number/, 'XYZ with a bad coordinate rejected'); throws(() => A.parseXYZ(''), /empty/, 'empty XYZ rejected');
  const pos5 = A.parsePOSCAR('ZnO wurtzite\n 3.2495\n 1.0 0.0 0.0\n -0.5 0.8660254038 0.0\n 0.0 0.0 1.6023696\n Zn_pv O\n 2 2\nSelective dynamics\nDirect\n 0.3333333 0.6666667 0.0 T T F\n 0.6666667 0.3333333 0.5 T T T\n 0.3333333 0.6666667 0.3819 T T T\n 0.6666667 0.3333333 0.8819 T T T\n');
  ok(pos5.atoms.length === 4 && pos5.atoms[0].el === 'Zn' && pos5.selectiveDynamics && pos5.atoms[0].fixed[2] === true, 'POSCAR VASP 5 with POTCAR-style labels and selective dynamics'); near(A.analyse(pos5).bonds[0].min_A, 1.974, 0.003, 'POSCAR Zn–O bond');
  const pos4 = A.parsePOSCAR('Si bulk\n-40.0475\n 0 .5 .5\n .5 0 .5\n .5 .5 0\n 2\nCartesian\n 0 0 0\n 0.25 0.25 0.25\n'); near(Math.hypot(pos4.atoms[1].x, pos4.atoms[1].y, pos4.atoms[1].z), 2.352, 0.05, 'POSCAR VASP 4, negative (volume) scale, Cartesian'); ok(pos4.atoms[0].el === 'Si', 'POSCAR VASP 4 species from the comment line');
  throws(() => A.parsePOSCAR('x\n1.0\n1 0 0\n0 1 0\n0 0 1\nSi\n2\nDirect\n0 0 0\n'), /coordinate lines|missing/, 'truncated POSCAR rejected'); throws(() => A.parsePOSCAR('x\nabc\n1 0 0\n0 1 0\n0 0 1\nSi\n1\nDirect\n0 0 0\n'), /scaling/, 'POSCAR bad scale rejected');
  const xsf = A.parseXSF('# comment\nCRYSTAL\nPRIMVEC\n 4.0 0 0\n 0 4.0 0\n 0 0 4.0\nCONVVEC\n 4.0 0 0\n 0 4.0 0\n 0 0 4.0\nPRIMCOORD\n 2 1\n 55 0 0 0 0.0 0.0 0.1\n Cl 2 2 2\n'); ok(xsf.atoms[0].el === 'Cs' && xsf.atoms[1].el === 'Cl' && xsf.cell[2][2] === 4, 'XSF with forces column and mixed species notation'); throws(() => A.parseXSF('CRYSTAL\nPRIMVEC\n1 0 0\n0 1 0\n0 0 1\n'), /PRIMCOORD/, 'XSF without atoms rejected');
  const pdb = A.parsePDB('HEADER    TEST\nCRYST1   10.000   12.000   15.000  90.00  90.00 120.00 P 1           1\nATOM      1  N   ALA A   1      11.104   6.134  -6.504  1.00  0.00           N  \nATOM      2  CA  ALA A   1      11.639   6.071  -5.147  1.00  0.00           C  \nHETATM    3 ZN    ZN A 101       1.000   2.000   3.000  1.00  0.00          ZN  \nHETATM    4 CA    CA A 102       4.000   5.000   6.000  1.00  0.00\nATOM      5  CA  GLY A   2       7.000   8.000   9.000\nEND\n');
  ok(pdb.atoms.map(a => a.el).join() === 'N,C,Zn,Ca,C' && pdb.atoms[0].z === -6.504, 'PDB element columns and name-based fallback (CA carbon vs calcium)'); near(A.cellParams(pdb.cell).gamma, 120, 1e-6, 'PDB CRYST1'); throws(() => A.parsePDB('REMARK nothing\nEND\n'), /no ATOM/, 'PDB without atoms rejected');
  const mol2 = A.parseMOL2('@<TRIPOS>MOLECULE\nethanol\n 3 2 1\nSMALL\nUSER_CHARGES\n\n@<TRIPOS>ATOM\n 1 C1 0.000 0.000 0.000 C.3 1 ETH -0.1\n 2 O1 1.430 0.000 0.000 O.3 1 ETH -0.6\n 3 CL1 -1.0 1.0 0.0 Cl 1 ETH 0.0\n@<TRIPOS>BOND\n 1 1 2 1\n 2 1 3 ar\n@<TRIPOS>CRYSIN\n 10 10 10 90 90 90 1 1\n'); ok(mol2.atoms.map(a => a.el).join() === 'C,O,Cl' && mol2.bonds.length === 2 && mol2.bonds[1][2] === 1.5 && mol2.title === 'ethanol' && mol2.cell[0][0] === 10 && mol2.atoms[1].q === -0.6, 'MOL2 atoms, types, bonds, charges, cell'); throws(() => A.parseMOL2('@<TRIPOS>MOLECULE\nx\n'), /ATOM/, 'MOL2 without atoms rejected');
  const sdf = A.parseSDF('methanol\n  test\n\n  2  1  0  0  0  0  0  0  0  0999 V2000\n    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n    1.4300    0.0000    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n  1  2  1  0\nM  CHG  1   2  -1\nM  END\n> <NAME>\nmethanol-ish\n\n$$$$\nsecond\n\n\n  1  0  0  0  0  0  0  0  0  0999 V2000\n    0.0000    0.0000    0.0000 Ar  0  0  0  0  0  0  0  0  0  0  0  0\nM  END\n$$$$\n');
  ok(sdf.atoms.length === 2 && sdf.bonds[0].join() === '0,1,1' && sdf.nMolecules === 2 && sdf.properties.NAME === 'methanol-ish' && sdf.atoms[1].charge === -1 && sdf.title === 'methanol', 'SDF first record, bonds, charges, data items'); throws(() => A.parseMOL('a\nb\nc\n  0  0 V3000\n'), /V3000/, 'V3000 molfile rejected with advice'); throws(() => A.parseMOL('a\nb\nc\n  5  0  0  0  0  0  0  0  0  0999 V2000\n    0.0000    0.0000    0.0000 C   0\n'), /declares 5 atoms/, 'truncated molfile rejected');
  const lmp = A.parseLAMMPSData('LAMMPS data file via hand\r\n\r\n4 atoms\r\n2 atom types\r\n\r\n-1.0 3.0 xlo xhi\r\n0.0 4.0 ylo yhi\r\n0.0 4.0 zlo zhi\r\n0.5 0.0 0.0 xy xz yz\r\n\r\nMasses\r\n\r\n1 132.905\r\n2 35.453 # Cl\r\n\r\nAtoms # full\r\n\r\n3 1 2 -1.0 1.0 2.0 2.0 0 0 0\r\n1 1 1 1.0 -1.0 0.0 0.0 0 0 0\r\n2 1 2 -1.0 1.0 2.0 2.0 0 0 0\r\n4 2 1 1.0 -1.0 0.0 0.0 1 0 0\r\n\r\nVelocities\r\n\r\n1 0 0 0\r\n\r\nBonds\r\n\r\n1 1 1 2\r\n');
  ok(lmp.atoms.map(a => a.el).join() === 'Cs,Cl,Cl,Cs' && lmp.atoms[0].x === 0 && lmp.atoms[1].q === -1 && lmp.cell[1][0] === 0.5 && lmp.cell[0][0] === 4 && lmp.atoms[3].mol === 2 && lmp.bonds.length === 1, 'LAMMPS data: full style, masses → elements, id sort, tilt, origin shift, image flags, CRLF');
  ok(A.parseLAMMPSData('t\n1 atoms\n1 atom types\n0 5 xlo xhi\n0 5 ylo yhi\n0 5 zlo zhi\n\nAtoms\n\n1 1 0.5 1 2 3\n', { elements: ['Na'] }).atoms[0].q === 0.5, 'LAMMPS data: charge style inferred, element map');
  throws(() => A.parseLAMMPSData('t\n2 atoms\n1 atom types\n0 5 xlo xhi\n0 5 ylo yhi\n0 5 zlo zhi\n\nAtoms # atomic\n\n1 1 0 0 0\n'), /declares 2 atoms/, 'LAMMPS data atom-count mismatch rejected'); throws(() => A.parseLAMMPSData('t\n1 atoms\n\nAtoms\n\n1 1 0 0 0\n'), /box bounds/, 'LAMMPS data without box rejected');
  const dump = 'ITEM: TIMESTEP\n0\nITEM: NUMBER OF ATOMS\n2\nITEM: BOX BOUNDS pp pp pp\n0 4\n0 4\n0 4\nITEM: ATOMS id type xs ys zs\n2 2 0.5 0.5 0.5\n1 1 0 0 0\nITEM: TIMESTEP\n100\nITEM: NUMBER OF ATOMS\n2\nITEM: BOX BOUNDS xy xz yz pp pp pp\n0 4.5 0.5\n0 4 0\n0 4 0\nITEM: ATOMS id type x y z\n1 1 0.1 0 0\n2 2 2.1 2 2\n';
  const d1 = A.parseLAMMPSDump(dump, { elements: ['Cs', 'Cl'] }), dA = A.parseLAMMPSDump(dump, { frames: 'all', elements: { 1: 'Cs', 2: 'Cl' } });
  ok(d1.frames.length === 1 && d1.atoms[1].el === 'Cl' && d1.atoms[1].x === 2 && d1.atoms[0].el === 'Cs', 'LAMMPS dump first frame, scaled coordinates, id sort'); ok(dA.frames.length === 2 && dA.frames[1].t === 100 && dA.frames[1].cell[1][0] === 0.5 && dA.frames[1].cell[0][0] === 4 && dA.frames[1].atoms[1].x === 2.1, 'LAMMPS dump all frames, triclinic bounds');
  ok(A.parseLAMMPSDump(dump).warnings.some(w => /elements/.test(w)), 'LAMMPS dump without element names warns'); throws(() => A.parseLAMMPSDump('nothing here'), /no "ITEM: ATOMS"/, 'non-dump text rejected');
  const gro = A.parseGRO('Zn in water\n    4\n    1SOL     OW    1   0.126   1.624   1.679  0.1227 -0.0580  0.0434\n    1SOL    HW1    2   0.190   1.661   1.747\n    2ZN      ZN    3   1.000   1.000   1.000\n    3ALA     CA    4   0.500   0.500   0.500\n   1.86206   1.86206   1.86206\n');
  ok(gro.atoms.map(a => a.el).join() === 'O,H,Zn,C' && Math.abs(gro.atoms[0].x - 1.26) < 1e-9 && Math.abs(gro.cell[0][0] - 18.6206) < 1e-9, 'GRO nm → Å, element inference from atom/residue names'); throws(() => A.parseGRO('t\n5\n    1SOL     OW    1   0.1   0.1   0.1\n 1 1 1\n'), /declares 5 atoms/, 'truncated GRO rejected');
  const cml = A.parseCML('<?xml version="1.0"?>\n<cml xmlns="http://www.xml-cml.org/schema">\n <molecule id="m1" title="water">\n  <atomArray>\n   <atom id="a1" elementType="O" x3="0.0" y3="0.0" z3="0.1173"/>\n   <atom id="a2" elementType="H" x3="0.0" y3="0.7572" z3="-0.4692"/>\n   <atom id=\'a3\' elementType=\'H\' x3=\'0.0\' y3=\'-0.7572\' z3=\'-0.4692\'/>\n  </atomArray>\n  <bondArray><bond atomRefs2="a1 a2" order="S"/><bond atomRefs2="a1 a3" order="1"/></bondArray>\n </molecule>\n</cml>');
  ok(cml.atoms.length === 3 && cml.bonds.length === 2 && cml.title === 'water' && cml.atoms[2].y === -0.7572, 'CML molecule with bonds');
  const cmlX = A.parseCML('<molecule><crystal><scalar title="a" units="units:angstrom">4.123(1)</scalar><scalar title="b">4.123</scalar><scalar title="c">4.123</scalar><scalar title="alpha">90</scalar><scalar title="beta">90</scalar><scalar title="gamma">90</scalar></crystal><atomArray><atom id="a1" elementType="Cs" xFract="0" yFract="0" zFract="0"/><atom id="a2" elementType="Cl" xFract="0.5" yFract="0.5" zFract="0.5"/></atomArray></molecule>'); near(cmlX.atoms[1].x, 2.0615, 1e-6, 'CML crystal with fractional coordinates'); throws(() => A.parseCML('<html><body/></html>'), /CML/, 'non-CML XML rejected');
  const gjf = A.parseGaussianInput('%chk=water.chk\n%mem=2GB\n# B3LYP/6-31G(d) opt\n  freq\n\nwater optimisation\n\n0 1\nO   0.000000   0.000000   0.117300\nH(Iso=2) 0 0.000000   0.757200  -0.469200\n1,0.0,-0.7572,-0.4692\n\n');
  ok(gjf.atoms.map(a => a.el).join('') === 'OHH' && gjf.charge === 0 && gjf.multiplicity === 1 && /B3LYP.*freq/.test(gjf.route) && gjf.atoms[1].y === 0.7572 && gjf.title === 'water optimisation', 'Gaussian input: Link 0, multi-line route, isotope tag, freeze flag, comma-separated record');
  const gpbc = A.parseGaussianInput('# pbepbe/3-21g\n\ncscl\n\n0 1\nCs 0.0 0.0 0.0\nCl 2.0 2.0 2.0\nTv 4.0 0.0 0.0\nTv 0.0 4.0 0.0\nTv 0.0 0.0 4.0\n'); ok(gpbc.atoms.length === 2 && gpbc.cell[2][2] === 4, 'Gaussian periodic input (Tv)');
  throws(() => A.parseGaussianInput('# hf/sto-3g\n\nzmat\n\n0 1\nO\nH 1 0.96\nH 1 0.96 2 104.5\n'), /Z-matrix/, 'Gaussian Z-matrix rejected with advice'); throws(() => A.parseGaussianInput('hello\n'), /route/, 'non-Gaussian text rejected');
  throws(() => A.parseGaussianInput('# hf/sto-3g\n\nzmat\n\n0 1\nO 0.0 0.0 0.0\nH 1 0.96 2 104.5\n'), /Z-matrix/, 'Gaussian mixed Z-matrix record rejected');
  const qe = A.parseQEInput("&CONTROL\n  calculation = 'scf' ! comment\n/\n&SYSTEM\n  ibrav = 0, nat = 4, ntyp = 2,\n  ecutwfc = 40.0\n/\n&ELECTRONS\n/\nATOMIC_SPECIES\nZn 65.38 Zn.pbe.UPF\nO1 15.999 O.pbe.UPF\nCELL_PARAMETERS (angstrom)\n  3.2495 0.0 0.0\n -1.62475 2.814149 0.0\n  0.0 0.0 5.2069\nATOMIC_POSITIONS {crystal}\nZn 0.3333333 0.6666667 0.0\nZn 0.6666667 0.3333333 0.5 0 0 0\nO1 1/3 2/3 0.3819\nO1 0.6666667 0.3333333 0.8819\nK_POINTS automatic\n6 6 4 0 0 0\n");
  ok(qe.atoms.length === 4 && A.analyse(qe).formula === 'ZnO', 'QE input: ibrav 0, crystal positions, labelled species, fractions'); near(A.analyse(qe).bonds[0].min_A, 1.974, 0.003, 'QE Zn–O bond');
  const qe2 = A.parseQEInput('&system\n ibrav=2, celldm(1)=10.2631, nat=2, ntyp=1\n/\nATOMIC_SPECIES\n Si 28.086 Si.UPF\nATOMIC_POSITIONS alat\n Si 0.00 0.00 0.00\n Si 0.25 0.25 0.25\nK_POINTS gamma\n'); near(A.analyse(qe2).nearestNeighbour_A, 2.352, 0.003, 'QE input: ibrav 2 with celldm(1) in bohr, alat positions'); near(A.analyse(qe2).density_kg_m3, 2329, 10, 'QE fcc primitive cell density');
  throws(() => A.parseQEInput('&system\n ibrav=5, celldm(1)=10, nat=1\n/\nATOMIC_POSITIONS crystal\n Si 0 0 0\n'), /ibrav = 5/, 'unsupported ibrav rejected with advice'); throws(() => A.parseQEInput('&system\n ibrav=0\n/\n'), /ATOMIC_POSITIONS/, 'QE input without positions rejected');
  throws(() => A.parseQEInput('&system\n ibrav=0, nat=3\n/\nCELL_PARAMETERS angstrom\n4 0 0\n0 4 0\n0 0 4\nATOMIC_POSITIONS angstrom\n Si 0 0 0\n'), /nat = 3/, 'QE nat mismatch rejected');
  const cst = A.parseCASTEPCell('# CASTEP cell\n%BLOCK LATTICE_ABC\nang\n 4.123 4.123 4.123\n 90 90 90\n%ENDBLOCK LATTICE_ABC\n\n%block positions_frac\nCs 0.0 0.0 0.0 ! comment\nCl:tag 0.5 0.5 0.5 SPIN=0.0\n%endblock positions_frac\nkpoint_mp_grid 4 4 4\n'); ok(cst.atoms.length === 2 && cst.atoms[1].el === 'Cl' && Math.abs(cst.atoms[1].z - 2.0615) < 1e-6, 'CASTEP cell: LATTICE_ABC, fractional positions, tags, comments');
  const cst2 = A.parseCASTEPCell('%BLOCK LATTICE_CART\nbohr\n 10 0 0\n 0 10 0\n 0 0 10\n%ENDBLOCK LATTICE_CART\n%BLOCK POSITIONS_ABS\nnm\nAr 0.1 0.2 0.3\n%ENDBLOCK POSITIONS_ABS\n'); near(cst2.cell[0][0], 5.29177, 1e-4, 'CASTEP cell: bohr lattice'); near(cst2.atoms[0].z, 3, 1e-9, 'CASTEP cell: nm absolute positions'); throws(() => A.parseCASTEPCell('%BLOCK POSITIONS_FRAC\nSi 0 0 0\n%ENDBLOCK POSITIONS_FRAC\n'), /LATTICE/, 'CASTEP cell without lattice rejected');
  for (const p of ['parseCIF', 'parseXYZ', 'parsePOSCAR', 'parseXSF', 'parsePDB', 'parseMOL', 'parseSDF', 'parseMOL2', 'parseLAMMPSData', 'parseGRO', 'parseCML', 'parseGaussianInput', 'parseQEInput', 'parseCASTEPCell', 'parseLAMMPSDump'])
    for (const junk of ['', 'garbage \u0000\u0001 text\nwith two lines\n', '{"json": true}']) { let good = false, msg = ''; try { A[p](junk); msg = 'returned a structure'; } catch (e) { good = e instanceof Error && e.message.length > 15 && !/undefined|is not a function|Cannot read/.test(e.message); msg = e.message; } ok(good, `${p} rejects junk input with a clear message`, msg); }

  // ── determinism, purity
  ok(JSON.stringify(A.nanoparticle('Si', { diameter_nm: 2 })) === JSON.stringify(A.nanoparticle('Si', { diameter_nm: 2 })), 'builders deterministic');
  const before = JSON.stringify(A.CRYSTALS.ZnO_wurtzite); A.slab('ZnO_wurtzite', {}); A.powderXRD('ZnO_wurtzite'); A.analyse(zno); ok(JSON.stringify(A.CRYSTALS.ZnO_wurtzite) === before, 'library data not mutated');

  log.unshift(`atomic geometry engine: ${pass} checks passed, ${fail} failed`);
  return { pass, fail, log };
}
