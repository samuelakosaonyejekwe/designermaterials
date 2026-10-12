// Suite 28 — Economics, Sustainability, Data and Lifecycle.
// Stoichiometric mass balance and techno-economic analysis, discounted cash flow with Monte Carlo and tornado sensitivity,
// screening cradle-to-gate life-cycle assessment and circularity, learning-curve and scale economics, and the digital
// material passport with end-of-life scenario comparison.
import * as N from '../core/num.js';
import * as chem from '../core/chem.js';

const F = N.fmt, OZ = 32.1507;   // troy ounces per kilogram
const kp = (key, label, value, unit, explain, good = null) => ({ key, label, value, unit, explain, good });
const pct = v => F(100 * v, 3) + ' %';
const money = v => (Math.abs(v) >= 1e6 ? F(v / 1e6, 3) + ' M' : Math.abs(v) >= 1e4 ? F(v / 1e3, 3) + ' k' : F(v, 3));

// ---------- production routes ----------
// re/pr: balanced with chem.js; buy: purchased form of each reactant (hydrate etc.); metal: element indexed to a market price
const ROUTES = {
  zno: { label: 'ZnO nanorods — hydrothermal, zinc nitrate + NaOH', re: [{ f: 'Zn(NO3)2', phase: 'aq', buy: 'Zn(NO3)2.6H2O' }, { f: 'NaOH', phase: 'aq' }], pr: [{ f: 'ZnO', phase: 's' }, { f: 'NaNO3', phase: 'aq' }, { f: 'H2O', phase: 'l' }], metal: 'Zn', name: 'Zinc oxide nanorods' },
  cuo: { label: 'CuO nanoparticles — precipitation, copper nitrate + NaOH', re: [{ f: 'Cu(NO3)2', phase: 'aq' }, { f: 'NaOH', phase: 'aq' }], pr: [{ f: 'CuO', phase: 's' }, { f: 'NaNO3', phase: 'aq' }, { f: 'H2O', phase: 'l' }], metal: 'Cu', name: 'Copper(II) oxide nanoparticles' },
  ag: { label: 'Ag nanoparticle colloid — AgNO₃ + NaBH₄', re: [{ f: 'AgNO3', phase: 'aq' }, { f: 'NaBH4', phase: 'aq' }, { f: 'H2O', phase: 'l' }], pr: [{ f: 'Ag', phase: 's' }, { f: 'H2', phase: 'g' }, { f: 'B(OH)3', phase: 'aq' }, { f: 'NaNO3', phase: 'aq' }], metal: 'Ag', name: 'Silver nanoparticles' },
  au: { label: 'Au nanoparticle colloid — HAuCl₄ + sodium citrate', re: [{ f: 'HAuCl4', phase: 'aq' }, { f: 'Na3C6H5O7', phase: 'aq' }], pr: [{ f: 'Au', phase: 's' }, { f: 'Na2C5H4O5', phase: 'aq' }, { f: 'CO2', phase: 'g' }, { f: 'NaCl', phase: 'aq' }, { f: 'HCl', phase: 'aq' }], metal: 'Au', name: 'Gold nanoparticles' },
  tio2: { label: 'TiO₂ nanoparticles / metasurface layer — alkoxide sol–gel', re: [{ f: 'Ti(OC3H7)4', phase: 'l' }, { f: 'H2O', phase: 'l' }], pr: [{ f: 'TiO2', phase: 's' }, { f: 'C3H8O', phase: 'l' }], metal: 'Ti', name: 'Titanium dioxide (sol–gel)' },
  resin: { label: 'Photopolymer lattice — additive manufacturing (no reaction stoichiometry)', re: [{ f: 'C5H8O2', phase: 'l' }], pr: [{ f: 'C5H8O2', phase: 's' }], metal: null, physical: true, name: 'Architected photopolymer lattice' },
};
const SOLVP = { water: { label: 'Water', rho: 1000, cp: 4180, hvap: 2.257e6, f: 'H2O' }, ethanol: { label: 'Ethanol', rho: 789, cp: 2440, hvap: 8.41e5, f: 'C2H6O' }, isopropanol: { label: 'Isopropanol', rho: 786, cp: 2600, hvap: 6.64e5, f: 'C3H8O' } };
// Elements that have appeared on EU critical-raw-material lists (screening only: verify against the current list)
const CRM = { Sb: 'antimony', Be: 'beryllium', Bi: 'bismuth', B: 'borates', Co: 'cobalt', Ga: 'gallium', Ge: 'germanium', Hf: 'hafnium', In: 'indium', Li: 'lithium', Mg: 'magnesium', Nb: 'niobium', Ta: 'tantalum', W: 'tungsten', V: 'vanadium', Sr: 'strontium', Pt: 'platinum-group metals', Pd: 'platinum-group metals', Rh: 'platinum-group metals', Ru: 'platinum-group metals', Ir: 'platinum-group metals', Os: 'platinum-group metals', Sc: 'scandium', Y: 'rare earths', La: 'rare earths', Ce: 'rare earths', Pr: 'rare earths', Nd: 'rare earths', Sm: 'rare earths', Eu: 'rare earths', Gd: 'rare earths', Tb: 'rare earths', Dy: 'rare earths', Ho: 'rare earths', Er: 'rare earths', Tm: 'rare earths', Yb: 'rare earths', Lu: 'rare earths', Ti: 'titanium metal (form-dependent)', Si: 'silicon metal (form-dependent)', Cu: 'copper (strategic, recent lists)', Ni: 'nickel (battery grade, recent lists)', P: 'phosphorus / phosphate rock' };
const LANG = { solid: 3.10, mixed: 3.63, fluid: 4.74 };
const EQUIP = [['Reactor and agitation', 0.30], ['Heating and utilities', 0.15], ['Solid–liquid separation and washing', 0.25], ['Drying and finishing', 0.15], ['Solvent recovery', 0.15]];

const SH = m => ({ mode: Array.isArray(m) ? m : [m] }), MB = ['tea', 'lca'];
const inputs = [
  { key: 'mode', label: 'Analysis', type: 'select', value: 'tea', group: 'Model', options: [
    { value: 'tea', label: 'Techno-economic analysis (mass balance, OPEX, CAPEX, unit cost)' }, { value: 'dcf', label: 'Discounted cash flow, Monte Carlo and tornado' }, { value: 'lca', label: 'Life-cycle screening, scenarios and circularity' },
    { value: 'learning', label: 'Learning curve, scale economics, make-or-buy' }, { value: 'passport', label: 'Digital material passport and end-of-life' }],
    help: 'Selects the economic or sustainability analysis. Cost and emission factors are user-adjustable screening values and must be replaced by quotations and qualified inventory data for decisions.' },
  { key: 'route', label: 'Production route', type: 'select', value: 'zno', group: 'Process', showIf: SH(['tea', 'lca', 'passport']), options: Object.entries(ROUTES).map(([value, r]) => ({ value, label: r.label })), help: 'Chemistry of the production process; the stoichiometry sets precursor demand, by-products and atom economy. The photopolymer route has no reaction stoichiometry and uses feedstock utilisation only.' },
  { key: 'cap', label: 'Annual production', unit: 'kg/yr', value: 20000, min: 1, max: 1e8, log: true, group: 'Process', showIf: SH(['tea', 'dcf']), help: 'Mass of saleable product per year at full capacity.' },
  { key: 'yieldFrac', label: 'Overall product yield', value: 0.9, min: 0.05, max: 1, step: 0.01, bus: 'fab.yield', group: 'Process', showIf: SH(MB), help: 'Saleable product divided by the stoichiometric maximum from the limiting precursor (conversion × isolation × quality yield).' },
  { key: 'excess', label: 'Excess of second reagent', unit: '%', value: 10, min: 0, max: 500, group: 'Process', showIf: SH(MB), help: 'Stoichiometric excess of the base or reducing agent over the balanced-equation requirement.' },
  { key: 'cPrec', label: 'Precursor concentration', unit: 'mol/m³', value: 200, min: 0.1, max: 5000, log: true, bus: 'proc.c_precursor', group: 'Process', showIf: SH(MB), help: 'Concentration of the limiting precursor in the reaction liquid (200 mol/m³ = 0.2 mol/L). Dilute recipes need large solvent volumes per kilogram of product.' },
  { key: 'washVol', label: 'Wash volume / reaction volume', value: 2, min: 0, max: 20, step: 0.5, group: 'Process', showIf: SH(MB), help: 'Solvent used for washing the product, as a multiple of the reaction volume.' },
  { key: 'solvType', label: 'Solvent', type: 'select', value: 'water', group: 'Process', showIf: SH(MB), options: Object.entries(SOLVP).map(([value, s]) => ({ value, label: s.label })), help: 'Process solvent; it sets density, heat capacity and the energy of recovery by distillation.' },
  { key: 'solvRec', label: 'Solvent recovery', unit: '%', value: 50, min: 0, max: 99, group: 'Process', showIf: SH(MB), help: 'Fraction of the used solvent recovered and reused; the remainder is make-up and becomes effluent.' },
  { key: 'T', label: 'Process temperature', unit: '°C', value: 90, min: 20, max: 300, bus: 'proc.T_C', group: 'Energy', showIf: SH(MB), help: 'Temperature to which the reaction liquid is heated from 20 °C.' },
  { key: 'heatLoss', label: 'Heat losses and hold', unit: '%', value: 25, min: 0, max: 300, group: 'Energy', showIf: SH(MB), help: 'Additional heat demand for losses during heat-up and hold, as a percentage of the sensible heat.' },
  { key: 'heatRec', label: 'Heat recovery', unit: '%', value: 0, min: 0, max: 90, group: 'Energy', showIf: SH(MB), help: 'Fraction of the sensible heat recovered from the discharged liquid to pre-heat the next feed.' },
  { key: 'eExtra', label: 'Other electricity (separation, drying, drives)', unit: 'kWh/kg', value: 3, min: 0, max: 5000, group: 'Energy', showIf: SH(MB), help: 'Electricity for stirring, pumping, centrifugation, drying or printing per kilogram of product.' },
  // prices
  { key: 'useIndex', label: 'Index precursor price to metal price', type: 'bool', value: true, group: 'Prices', showIf: SH('tea'), help: 'When on, the precursor price is metal content × metal market price × premium, so it follows the live metal price. When off, or for routes without an indexed metal, the entered precursor price is used.' },
  { key: 'premium', label: 'Precursor premium over metal value', value: 4, min: 1, max: 50, step: 0.1, group: 'Prices', showIf: SH('tea'), help: 'Ratio of the precursor price to the market value of the metal it contains (conversion, purification and distribution): about 3–6 for zinc and copper salts, 1.1–1.5 for gold and silver salts.' },
  { key: 'pPrec', label: 'Precursor price (not indexed)', unit: 'USD/kg', value: 2.5, min: 0.01, max: 1e6, log: true, group: 'Prices', showIf: SH('tea'), help: 'Price of the metal precursor or feedstock as purchased, used when indexing is off or no metal index applies (titanium alkoxide, photopolymer resin).' },
  { key: 'pReag', label: 'Second-reagent price', unit: 'USD/kg', value: 0.5, min: 0, max: 1e5, group: 'Prices', showIf: SH('tea'), help: 'Price of the base, reducing agent or other co-reagent.' },
  { key: 'pSolv', label: 'Solvent price', unit: 'USD/m³', value: 2, min: 0, max: 1e5, group: 'Prices', showIf: SH('tea'), help: 'Price of fresh solvent: deionised water a few USD/m³, technical ethanol or isopropanol of the order of 1000 USD/m³.' },
  { key: 'zinc', label: 'Zinc price', unit: 'USD/t', value: 2800, min: 500, max: 20000, live: 'price.zinc_usd_t', group: 'Market data', showIf: SH(['tea', 'passport']), help: 'Market price of zinc metal; the default is indicative, use the live feed for a current value.' },
  { key: 'copper', label: 'Copper price', unit: 'USD/t', value: 9800, min: 1000, max: 50000, live: 'price.copper_usd_t', group: 'Market data', showIf: SH(['tea', 'passport']), help: 'Market price of copper metal; indicative default, use the live feed.' },
  { key: 'silver', label: 'Silver price', unit: 'USD/oz', value: 38, min: 5, max: 500, live: 'price.silver_usd_oz', group: 'Market data', showIf: SH(['tea', 'passport']), help: 'Market price of silver per troy ounce; indicative default, use the live feed.' },
  { key: 'gold', label: 'Gold price', unit: 'USD/oz', value: 3300, min: 500, max: 20000, live: 'price.gold_usd_oz', group: 'Market data', showIf: SH(['tea', 'passport']), help: 'Market price of gold per troy ounce; indicative default, use the live feed.' },
  { key: 'fx', label: 'Exchange rate', unit: 'USD/EUR', value: 1.10, min: 0.5, max: 2, step: 0.01, live: 'fx.USD_per_EUR', group: 'Market data', showIf: SH(['tea', 'dcf']), help: 'US dollars per euro, used to convert the electricity price and to report results in both currencies.' },
  { key: 'pElec', label: 'Electricity price', unit: 'EUR/kWh', value: 0.15, min: 0.01, max: 2, step: 0.01, group: 'Prices', showIf: SH('tea'), help: 'Industrial electricity price including grid charges.' },
  { key: 'nOp', label: 'Operators per shift', value: 1, min: 0, max: 200, step: 1, group: 'Fixed costs', showIf: SH('tea'), help: 'Number of operators present during operating hours.' },
  { key: 'wage', label: 'Loaded labour rate', unit: 'USD/h', value: 35, min: 1, max: 300, group: 'Fixed costs', showIf: SH('tea'), help: 'Hourly labour cost including social charges.' },
  { key: 'opHours', label: 'Operating hours per year', unit: 'h/yr', value: 6000, min: 200, max: 8760, group: 'Fixed costs', showIf: SH('tea'), help: 'Manned operating time per year.' },
  { key: 'overhead', label: 'Overhead on labour', unit: '%', value: 30, min: 0, max: 300, group: 'Fixed costs', showIf: SH('tea'), help: 'Supervision, laboratory, administration and insurance, as a percentage of direct labour.' },
  { key: 'maint', label: 'Maintenance', unit: '% of CAPEX/yr', value: 4, min: 0, max: 20, group: 'Fixed costs', showIf: SH('tea'), help: 'Annual maintenance and repair as a percentage of fixed capital investment.' },
  { key: 'pWaste', label: 'Solid/salt waste treatment', unit: 'USD/kg', value: 0.15, min: 0, max: 100, group: 'Prices', showIf: SH('tea'), help: 'Cost of treating or disposing of by-product salts, unreacted reagents and unrecovered organic solvent.' },
  { key: 'pWW', label: 'Wastewater treatment', unit: 'USD/m³', value: 3, min: 0, max: 500, group: 'Prices', showIf: SH('tea'), help: 'Cost of treating aqueous effluent.' },
  { key: 'eqRef', label: 'Purchased equipment cost at reference capacity', unit: 'USD', value: 400000, min: 1000, max: 1e9, log: true, group: 'Capital', showIf: SH('tea'), help: 'Total purchased cost of the main equipment for a plant of the reference capacity, from quotations or cost charts.' },
  { key: 'capRef', label: 'Reference capacity', unit: 'kg/yr', value: 20000, min: 1, max: 1e8, log: true, group: 'Capital', showIf: SH('tea'), help: 'Capacity to which the reference equipment cost refers.' },
  { key: 'scaleExp', label: 'Capacity exponent', value: 0.6, min: 0.3, max: 1, step: 0.05, group: 'Capital', showIf: SH(['tea', 'learning']), help: 'Exponent of the cost–capacity law C = C_ref·(S/S_ref)ⁿ; 0.6 is the classical six-tenths rule, numbering-up of identical units approaches 1.' },
  { key: 'lang', label: 'Lang factor (plant type)', type: 'select', value: 'mixed', group: 'Capital', showIf: SH('tea'), options: [{ value: 'solid', label: 'Solid processing (3.10)' }, { value: 'mixed', label: 'Solid–fluid processing (3.63)' }, { value: 'fluid', label: 'Fluid processing (4.74)' }], help: 'Multiplier from purchased equipment cost to fixed capital investment (installation, piping, instrumentation, buildings, engineering).' },
  { key: 'life', label: 'Project life', unit: 'yr', value: 10, min: 2, max: 40, step: 1, group: 'Finance', showIf: SH(['tea', 'dcf']), help: 'Economic life over which the investment is recovered.' },
  { key: 'rate', label: 'Discount rate', unit: '%', value: 10, min: 0, max: 40, group: 'Finance', showIf: SH(['tea', 'dcf']), help: 'Required return on capital used for discounting and for the capital recovery factor.' },
  { key: 'price', label: 'Selling price', unit: 'USD/kg', value: 60, min: 0.01, max: 1e7, log: true, group: 'Finance', showIf: SH(['tea', 'dcf']), help: 'Net selling price of the product.' },
  { key: 'loading', label: 'Areal loading in the application', unit: 'g/m²', value: 2, min: 1e-4, max: 1e5, log: true, group: 'Application', showIf: SH(MB), help: 'Mass of the material per square metre of coated or patterned surface, for cost and impact per square metre.' },
  { key: 'devArea', label: 'Device area', unit: 'cm²', value: 1, min: 1e-4, max: 1e5, log: true, group: 'Application', showIf: SH(MB), help: 'Functional area of one device, for cost and impact per device.' },
  { key: 'integ', label: 'Integration cost per device', unit: 'USD', value: 0.05, min: 0, max: 1e5, group: 'Application', showIf: SH('tea'), help: 'Assembly, patterning and packaging cost added per device, excluding the material itself.' },
  // DCF
  { key: 'capex', label: 'Capital investment', unit: 'USD', value: 1.45e6, min: 1000, max: 1e10, log: true, bus: 'econ.capex', group: 'Cash flow', showIf: SH('dcf'), help: 'Fixed capital investment spent in year 0 (from the techno-economic analysis).' },
  { key: 'opex', label: 'Operating cost at full capacity', unit: 'USD/yr', value: 6.0e5, min: 0, max: 1e10, bus: 'econ.opex', group: 'Cash flow', showIf: SH('dcf'), help: 'Annual operating expenditure at full capacity, excluding depreciation and capital charges.' },
  { key: 'varFrac', label: 'Variable share of operating cost', unit: '%', value: 45, min: 0, max: 100, group: 'Cash flow', showIf: SH('dcf'), help: 'Share of the operating cost that scales with output (materials, energy, waste); the rest is fixed (labour, maintenance).' },
  { key: 'ramp', label: 'Ramp-up period', unit: 'yr', value: 2, min: 0, max: 10, step: 1, group: 'Cash flow', showIf: SH('dcf'), help: 'Years to reach full output; sales grow linearly during this period.' },
  { key: 'tax', label: 'Income tax rate', unit: '%', value: 25, min: 0, max: 60, group: 'Cash flow', showIf: SH('dcf'), help: 'Tax on positive taxable income (revenue − operating cost − depreciation). Losses are not carried forward in this model.' },
  { key: 'depYears', label: 'Depreciation period', unit: 'yr', value: 7, min: 1, max: 40, step: 1, group: 'Cash flow', showIf: SH('dcf'), help: 'Straight-line depreciation period of the capital investment for tax purposes.' },
  { key: 'wc', label: 'Working capital', unit: '% of revenue', value: 10, min: 0, max: 50, group: 'Cash flow', showIf: SH('dcf'), help: 'Inventory and receivables tied up, as a percentage of annual revenue; recovered at the end of the project.' },
  { key: 'cvPrice', label: 'Uncertainty of price', unit: '% CV', value: 15, min: 0, max: 80, group: 'Uncertainty', showIf: SH('dcf'), help: 'Coefficient of variation of the selling price in the Monte Carlo analysis (lognormal, constant over the project).' },
  { key: 'cvVol', label: 'Uncertainty of sales volume', unit: '% CV', value: 15, min: 0, max: 80, group: 'Uncertainty', showIf: SH('dcf'), help: 'Coefficient of variation of the sales volume.' },
  { key: 'cvOpex', label: 'Uncertainty of operating cost', unit: '% CV', value: 10, min: 0, max: 80, group: 'Uncertainty', showIf: SH('dcf'), help: 'Coefficient of variation of the operating cost.' },
  { key: 'cvCapex', label: 'Uncertainty of capital cost', unit: '% CV', value: 20, min: 0, max: 80, group: 'Uncertainty', showIf: SH('dcf'), help: 'Coefficient of variation of the capital investment; factored estimates of this kind are typically accurate to ±30 %.' },
  { key: 'nMc', label: 'Monte Carlo samples', value: 3000, min: 200, max: 30000, step: 100, sens: false, group: 'Uncertainty', showIf: SH('dcf'), help: 'Number of random scenarios for the NPV distribution.' },
  { key: 'swing', label: 'Tornado swing', unit: '%', value: 20, min: 1, max: 80, group: 'Uncertainty', showIf: SH('dcf'), help: 'One-at-a-time variation applied to each input for the tornado sensitivity table.' },
  { key: 'seed', label: 'Random seed', value: 1, min: 1, max: 99999, step: 1, sens: false, group: 'Uncertainty', showIf: SH('dcf'), help: 'Seed of the pseudo-random generator; the same seed reproduces the same Monte Carlo result.' },
  // LCA
  { key: 'fu', label: 'Functional unit', type: 'select', value: 'kg', group: 'Life cycle', showIf: SH('lca'), options: [{ value: 'kg', label: '1 kg of product at the factory gate' }, { value: 'm2', label: '1 m² of coated or patterned surface' }, { value: 'device', label: '1 device' }], help: 'Reference quantity to which all impacts are normalised; it must reflect the function delivered for comparisons to be meaningful.' },
  { key: 'efPrec', label: 'Emission factor: precursor', unit: 'kg CO₂e/kg', value: 3, min: 0, max: 50000, group: 'Screening factors', showIf: SH('lca'), help: 'Cradle-to-gate greenhouse-gas emissions of producing the precursor as purchased. User-adjustable screening factor: base metal salts are of the order of 1–5, silver salts hundreds, gold salts tens of thousands of kg CO₂e/kg.' },
  { key: 'efReag', label: 'Emission factor: second reagent', unit: 'kg CO₂e/kg', value: 1.2, min: 0, max: 1000, group: 'Screening factors', showIf: SH('lca'), help: 'Cradle-to-gate emissions of the co-reagent (sodium hydroxide is of the order of 1 kg CO₂e/kg). User-adjustable screening factor.' },
  { key: 'efSolv', label: 'Emission factor: solvent', unit: 'kg CO₂e/kg', value: 0.0005, min: 0, max: 100, group: 'Screening factors', showIf: SH('lca'), help: 'Cradle-to-gate emissions of fresh solvent: deionised water below 0.001, fossil-based alcohols of the order of 1–2 kg CO₂e/kg. User-adjustable screening factor.' },
  { key: 'grid', label: 'Grid carbon intensity', unit: 'g CO₂/kWh', value: 350, min: 0, max: 1200, live: 'grid.co2_g_kWh', group: 'Screening factors', showIf: SH(['lca', 'passport']), help: 'Carbon intensity of the electricity used; take the live value for the plant location.' },
  { key: 'gridRen', label: 'Renewable-supply carbon intensity', unit: 'g CO₂/kWh', value: 30, min: 0, max: 300, group: 'Screening factors', showIf: SH('lca'), help: 'Carbon intensity assumed for the renewable-electricity scenario (life-cycle value of wind and solar supply).' },
  { key: 'efWaste', label: 'Emission factor: waste treatment', unit: 'kg CO₂e/kg', value: 0.1, min: 0, max: 20, group: 'Screening factors', showIf: SH('lca'), help: 'Emissions of treating by-product salts and residues. User-adjustable screening factor.' },
  { key: 'transKm', label: 'Transport distance of inputs', unit: 'km', value: 500, min: 0, max: 30000, group: 'Screening factors', showIf: SH('lca'), help: 'Average road distance over which precursors and reagents are delivered.' },
  { key: 'efTrans', label: 'Emission factor: road freight', unit: 'kg CO₂e/t·km', value: 0.1, min: 0, max: 2, step: 0.01, group: 'Screening factors', showIf: SH('lca'), help: 'Greenhouse-gas emissions per tonne-kilometre of road freight. User-adjustable screening factor.' },
  { key: 'cedPrec', label: 'Embodied energy: precursor', unit: 'MJ/kg', value: 40, min: 0, max: 1e6, group: 'Screening factors', showIf: SH('lca'), help: 'Cumulative primary energy demand of producing the precursor. User-adjustable screening factor.' },
  { key: 'cedReag', label: 'Embodied energy: second reagent', unit: 'MJ/kg', value: 20, min: 0, max: 1e5, group: 'Screening factors', showIf: SH('lca'), help: 'Cumulative primary energy demand of the co-reagent. User-adjustable screening factor.' },
  { key: 'pef', label: 'Primary energy factor of electricity', unit: 'MJ/kWh', value: 8, min: 3.6, max: 15, step: 0.1, group: 'Screening factors', showIf: SH('lca'), help: 'Primary energy consumed per kilowatt-hour of electricity delivered (3.6 = ideal conversion; fossil-dominated grids 9–11).' },
  { key: 'coolWater', label: 'Cooling and utility water', unit: 'L/kg', value: 50, min: 0, max: 1e5, group: 'Screening factors', showIf: SH('lca'), help: 'Water consumed for cooling and utilities per kilogram of product, in addition to process water.' },
  { key: 'eSave', label: 'Energy saved in use', unit: 'kWh/yr per kg', value: 50, min: 0, max: 1e6, group: 'Application', showIf: SH('lca'), help: 'Electricity saved per year in the application by each kilogram of deployed material (for example a radiative-cooling or low-emissivity metasurface); used for the carbon payback time.' },
  { key: 'recIn', label: 'Recycled content of feedstock', unit: '%', value: 10, min: 0, max: 100, group: 'Circularity', showIf: SH('lca'), help: 'Share of the material input that comes from recycled sources.' },
  { key: 'recOut', label: 'Collected for recycling at end of life', unit: '%', value: 30, min: 0, max: 100, group: 'Circularity', showIf: SH(['lca', 'passport']), help: 'Share of the product collected and sent to recycling after use.' },
  { key: 'recEff', label: 'Recycling process efficiency', unit: '%', value: 80, min: 1, max: 100, group: 'Circularity', showIf: SH(['lca', 'passport']), help: 'Fraction of the collected material actually recovered by the recycling process.' },
  { key: 'lifeRatio', label: 'Product lifetime / industry average', value: 1, min: 0.1, max: 10, step: 0.1, group: 'Circularity', showIf: SH('lca'), help: 'Utility factor of the circularity indicator: lifetime (or intensity of use) relative to an average product of the same type.' },
  { key: 'elements', label: 'Additional formula for critical-material screening', type: 'text', value: '', group: 'Circularity', showIf: SH(['lca', 'passport']), help: 'Optional chemical formula of dopants, electrodes or substrates in the product (for example In2O3 or GaN) to include in the critical-raw-material screening.' },
  // learning
  { key: 'c0', label: 'Present unit cost', unit: 'USD/kg', value: 45, min: 0.01, max: 1e7, log: true, bus: 'fab.cost_per_kg', group: 'Learning curve', showIf: SH('learning'), help: 'Unit production cost at the present cumulative volume.' },
  { key: 'q0', label: 'Cumulative production to date', unit: 'kg', value: 5000, min: 1, max: 1e10, log: true, group: 'Learning curve', showIf: SH('learning'), help: 'Total quantity produced so far, the reference point of the experience curve.' },
  { key: 'lr', label: 'Learning rate', unit: '% per doubling', value: 15, min: 0, max: 50, group: 'Learning curve', showIf: SH('learning'), help: 'Cost reduction for each doubling of cumulative production. Process industries commonly show 10–25 %; the value must come from comparable technologies.' },
  { key: 'qT', label: 'Target cumulative production', unit: 'kg', value: 1e6, min: 1, max: 1e12, log: true, group: 'Learning curve', showIf: SH('learning'), help: 'Cumulative volume at which the future cost is evaluated.' },
  { key: 'cTarget', label: 'Target unit cost', unit: 'USD/kg', value: 15, min: 0.001, max: 1e7, log: true, group: 'Learning curve', showIf: SH('learning'), help: 'Cost needed for the application to be competitive.' },
  { key: 'qYear', label: 'Present annual volume', unit: 'kg/yr', value: 2000, min: 0.1, max: 1e9, log: true, group: 'Learning curve', showIf: SH('learning'), help: 'Production in the current year.' },
  { key: 'growth', label: 'Annual volume growth', unit: '%/yr', value: 40, min: 0, max: 300, group: 'Learning curve', showIf: SH('learning'), help: 'Growth rate of annual production used to translate cumulative volume into calendar time.' },
  { key: 'sRatio', label: 'Plant scale-up factor', value: 10, min: 0.1, max: 1000, log: true, group: 'Scale economics', showIf: SH('learning'), help: 'Ratio of the new plant capacity to the present one.' },
  { key: 'capShare', label: 'Capital-related share of unit cost', unit: '%', value: 35, min: 0, max: 100, group: 'Scale economics', showIf: SH('learning'), help: 'Share of the present unit cost that scales with capital (depreciation, maintenance, largely also labour).' },
  { key: 'fixedMake', label: 'Make: annual fixed cost', unit: 'USD/yr', value: 250000, min: 0, max: 1e9, group: 'Make or buy', showIf: SH('learning'), help: 'Annual fixed cost of in-house production (capital charge, labour, maintenance).' },
  { key: 'varMake', label: 'Make: variable cost', unit: 'USD/kg', value: 20, min: 0, max: 1e6, group: 'Make or buy', showIf: SH('learning'), help: 'Variable cost of in-house production per kilogram (materials, energy, waste).' },
  { key: 'buyPrice', label: 'Buy: purchase price', unit: 'USD/kg', value: 55, min: 0.01, max: 1e7, log: true, group: 'Make or buy', showIf: SH('learning'), help: 'Delivered price of the material from an external supplier.' },
  { key: 'demandMB', label: 'Annual requirement', unit: 'kg/yr', value: 8000, min: 0.1, max: 1e9, log: true, group: 'Make or buy', showIf: SH('learning'), help: 'Quantity needed per year.' },
  // passport
  { key: 'pName', label: 'Product name', type: 'text', value: '', group: 'Identity', showIf: SH('passport'), bus: 'material.name', help: 'Trade or descriptive name of the material; if empty the name of the selected route is used.' },
  { key: 'pFormula', label: 'Chemical formula', type: 'text', value: '', group: 'Identity', showIf: SH('passport'), bus: 'material.formula', help: 'Chemical formula of the product; if empty the product of the selected route is used.' },
  { key: 'batch', label: 'Batch or lot identifier', type: 'text', value: '', group: 'Provenance', showIf: SH('passport'), help: 'Identifier that links the passport to production records.' },
  { key: 'producer', label: 'Producer', type: 'text', value: '', group: 'Provenance', showIf: SH('passport'), help: 'Legal entity and site that manufactured the batch.' },
  { key: 'issued', label: 'Date of issue', type: 'text', value: '', group: 'Provenance', showIf: SH('passport'), help: 'Date of issue of the passport (ISO format, e.g. 2026-10-12).' },
  { key: 'hazText', label: 'Hazard classification', type: 'text', value: '', group: 'Identity', showIf: SH('passport'), help: 'Hazard statements and pictograms from the safety data sheet of the material. Leave empty if not yet established; the passport then reports the field as missing.' },
  { key: 'eolRoute', label: 'Declared end-of-life route', type: 'select', value: 'recycling', group: 'End of life', showIf: SH('passport'), options: [{ value: 'recycling', label: 'Material recycling' }, { value: 'incineration', label: 'Incineration with energy recovery' }, { value: 'landfill', label: 'Landfill' }], help: 'Route recorded in the passport; all three are compared on cost and greenhouse-gas emissions.' },
  { key: 'co2Primary', label: 'Carbon footprint of primary production', unit: 'kg CO₂e/kg', value: 12, min: 0, max: 1e5, bus: 'econ.co2_per_kg', group: 'End of life', showIf: SH('passport'), help: 'Cradle-to-gate footprint of the product (from the life-cycle screening); recycling is credited with the primary production it avoids.' },
  { key: 'feeLand', label: 'Landfill gate fee', unit: 'USD/t', value: 120, min: 0, max: 5000, group: 'End of life', showIf: SH('passport'), help: 'Disposal charge for landfilling the end-of-life product.' },
  { key: 'feeInc', label: 'Incineration gate fee', unit: 'USD/t', value: 180, min: 0, max: 5000, group: 'End of life', showIf: SH('passport'), help: 'Charge for incineration with energy recovery.' },
  { key: 'costRec', label: 'Recycling process cost', unit: 'USD/kg', value: 2, min: 0, max: 1e5, group: 'End of life', showIf: SH('passport'), help: 'Collection, separation and reprocessing cost per kilogram of end-of-life product.' },
  { key: 'eRec', label: 'Recycling process energy', unit: 'kWh/kg', value: 3, min: 0, max: 1000, group: 'End of life', showIf: SH('passport'), help: 'Electricity used by the recycling process per kilogram treated.' },
  { key: 'valRec', label: 'Value of recovered material', unit: 'USD/kg', value: 0, min: 0, max: 1e6, group: 'End of life', showIf: SH('passport'), help: 'Market value of the recovered material; enter 0 to use the indexed metal value of the product where a metal price is available.' },
  { key: 'lhv', label: 'Heating value of the product', unit: 'MJ/kg', value: 0, min: 0, max: 50, group: 'End of life', showIf: SH('passport'), help: 'Lower heating value released on incineration: zero for oxides and metals, about 25–30 MJ/kg for acrylate photopolymers.' },
];

// ---------- shared: stoichiometric mass and energy inventory per kg of product ----------
function inventory(inp, over = {}) {
  const q = { ...inp, ...over }, rt = ROUTES[q.route] || ROUTES.zno, sv = SOLVP[q.solvType] || SOLVP.water, re = rt.re.map(s => s.f), pr = rt.pr.map(s => s.f), y = N.clamp(q.yieldFrac, 0.01, 1);
  const Mp = chem.molarMass(pr[0]) / 1e3, nP = 1 / Mp; let bal = null, rows = [], mPrec, mReag = 0, molPrec, byprod = 0, unreacted, atomEco = 1;
  if (rt.physical) { mPrec = 1 / y; molPrec = mPrec / Mp; unreacted = mPrec - 1; rows = [[chem.prettyFormula(re[0]), 'feedstock (resin)', 1, mPrec], [chem.prettyFormula(pr[0]) + ' (cured)', 'product', 1, 1]]; }
  else {
    bal = chem.balance(re, pr); const nuP = bal.products[0], x = nP / nuP;   // reaction extent per kg product (mol), before yield
    molPrec = x * bal.reactants[0] / y; mPrec = molPrec * chem.molarMass(rt.re[0].buy || re[0]) / 1e3; unreacted = (1 - y) * molPrec * chem.molarMass(re[0]) / 1e3;
    rows.push([chem.prettyFormula(rt.re[0].buy || re[0]), 'precursor', bal.reactants[0], mPrec]);
    rt.re.slice(1).forEach((s, i) => { const isW = s.f === 'H2O', mol = x * bal.reactants[i + 1] / y * (isW ? 1 : 1 + q.excess / 100), m = mol * chem.molarMass(s.f) / 1e3; if (!isW) { mReag += m; unreacted += (mol - x * bal.reactants[i + 1]) * chem.molarMass(s.f) / 1e3; } rows.push([chem.prettyFormula(s.f), isW ? 'reactant (from solvent)' : 'reagent', bal.reactants[i + 1], m]); });
    rows.push([chem.prettyFormula(pr[0]), 'product', nuP, 1]);
    rt.pr.slice(1).forEach((s, i) => { const m = x * bal.products[i + 1] * chem.molarMass(s.f) / 1e3; if (s.f !== 'H2O') byprod += m; rows.push([chem.prettyFormula(s.f), 'by-product', bal.products[i + 1], m]); });
    atomEco = nuP * chem.molarMass(pr[0]) / N.sum(re.map((f, i) => bal.reactants[i] * chem.molarMass(f)));
  }
  const Vrxn = rt.physical ? 0 : molPrec / q.cPrec, Vwash = rt.physical ? q.washVol * 2e-3 : q.washVol * Vrxn, Vtot = Vrxn + Vwash, rec = q.solvRec / 100, mSolv = Vtot * sv.rho, mMake = mSolv * (1 - rec), organic = q.solvType !== 'water';
  const dT = Math.max(q.T - 20, 0), eHeat = sv.rho * sv.cp * Vrxn * dT * (1 + q.heatLoss / 100) * (1 - q.heatRec / 100) / 3.6e6, eRecov = organic ? rec * mSolv * sv.hvap * 1.3 / 3.6e6 : rec * Vtot * 0.5, eTot = eHeat + eRecov + q.eExtra;
  const waste = byprod + unreacted + (organic ? mMake : 0), ww = organic ? 0 : Vtot * (1 - rec), eFactor = waste + (rt.physical ? 0 : 0);
  return { rt, sv, re, pr, bal, rows, Mp, mPrec, mReag, molPrec, byprod, unreacted, atomEco, Vrxn, Vtot, mSolv, mMake, eHeat, eRecov, eTot, waste, ww, eFactor, organic, y };
}
const metalPrice = (inp, el) => ({ Zn: inp.zinc / 1e3, Cu: inp.copper / 1e3, Ag: inp.silver * OZ, Au: inp.gold * OZ }[el] ?? null);   // USD per kg metal
function routeEquations(iv) {
  const rt = iv.rt; if (rt.physical) return [{ kind: 'formula', label: 'Feedstock', text: `n ${chem.prettyFormula(iv.re[0])} (liquid resin) → (${chem.prettyFormula(iv.re[0])})ₙ (cured network)`, note: 'Representative acrylate repeat unit; photopolymerisation has no by-products, so the mass balance is feedstock utilisation only.' }];
  const ion = chem.ionicEquations(rt.re, rt.pr), out = [{ kind: 'chemical', label: 'Balanced production reaction', text: ion.molecular }];
  if (ion.completeIonic !== ion.molecular) out.push({ kind: 'ionic', label: 'Complete ionic equation', text: ion.completeIonic }, { kind: 'ionic', label: 'Net ionic equation', text: ion.netIonic, note: ion.spectators.length ? 'Spectator ions: ' + [...new Set(ion.spectators)].join(', ') : undefined });
  return out;
}
const crf = (i, n) => i > 1e-9 ? i * Math.pow(1 + i, n) / (Math.pow(1 + i, n) - 1) : 1 / n;

// ---------- mode: techno-economic analysis ----------
function runTea(inp) {
  const w = [], iv = inventory(inp), rt = iv.rt, cap = inp.cap, wMetal = rt.metal && !rt.physical ? (chem.massFractions(rt.re[0].buy || iv.re[0])[rt.metal] || 0) : 0, pm = rt.metal ? metalPrice(inp, rt.metal) : null;
  const indexed = inp.useIndex && pm !== null && wMetal > 0, pPrec = indexed ? wMetal * pm * inp.premium : inp.pPrec;
  if (inp.useIndex && !indexed) w.push(`No market-price index is available for this route (${rt.metal || 'polymer feedstock'}): the entered precursor price of ${inp.pPrec} USD/kg is used.`);
  const cMat = iv.mPrec * pPrec + iv.mReag * inp.pReag, cSolv = iv.mMake / iv.sv.rho * inp.pSolv, pEl = inp.pElec * inp.fx, cEn = iv.eTot * pEl, cWaste = iv.waste * inp.pWaste + iv.ww * inp.pWW, varKg = cMat + cSolv + cEn + cWaste;
  const used = EQUIP.filter(e => e[0] !== 'Solvent recovery' || inp.solvRec > 0), sc = Math.pow(cap / inp.capRef, inp.scaleExp), eq = used.map(([nm, f]) => [nm, f * inp.eqRef, f * inp.eqRef * sc]), purch = N.sum(eq.map(e => e[2])), LF = LANG[inp.lang] || 3.63, capex = LF * purch;
  const labour = inp.nOp * inp.wage * inp.opHours * (1 + inp.overhead / 100), maint = inp.maint / 100 * capex, CRF = crf(inp.rate / 100, inp.life), capCh = CRF * capex, fixed = labour + maint, opex = varKg * cap + fixed;
  const unit = (opex + capCh) / cap, margin = (inp.price - unit) / inp.price, contrib = inp.price - varKg, be = contrib > 0 ? (fixed + capCh) / contrib : 0, perM2 = unit * inp.loading / 1e3, perDev = perM2 * inp.devArea / 1e4 + inp.integ;
  if (contrib <= 0) w.push(`The selling price (${inp.price} USD/kg) does not cover the variable cost (${F(varKg, 3)} USD/kg): there is no break-even volume; it is reported as 0.`);
  if (be > cap) w.push(`Break-even volume (${F(be, 4)} kg/yr) exceeds the plant capacity: the project loses money at full output with these inputs.`);
  if (cap / inp.capRef > 10 || cap / inp.capRef < 0.1) w.push('Capacity differs from the reference by more than a factor of ten: the cost–capacity exponent is being extrapolated beyond its usual range of validity.');
  const cats = [['Precursor', iv.mPrec * pPrec], ['Second reagent', iv.mReag * inp.pReag], ['Solvent make-up', cSolv], ['Energy', cEn], ['Waste and effluent', cWaste], ['Labour and overhead', labour / cap], ['Maintenance', maint / cap], ['Capital charge', capCh / cap]], top = cats.slice().sort((a, b) => b[1] - a[1])[0];
  const caps = N.logspace(Math.log10(cap / 20), Math.log10(cap * 20), 50), unitAt = c => { const cx = LF * N.sum(used.map(([, f]) => f * inp.eqRef)) * Math.pow(c / inp.capRef, inp.scaleExp); return varKg + (labour + inp.maint / 100 * cx + CRF * cx) / c; };
  const ms = N.linspace(0.5, 2, 31), unitMetal = k => unit + iv.mPrec * pPrec * (k - 1), eqs = routeEquations(iv);
  eqs.push({ kind: 'formula', label: 'Precursor demand and atom economy', text: 'm_prec = (ν_prec·M_prec)/(ν_prod·M_prod·Y);   atom economy = ν_prod·M_prod/Σ(νᵢ·Mᵢ)_reactants' },
    { kind: 'formula', label: 'Operating and capital cost', text: 'OPEX = Σ(mᵢ·pᵢ) + E·p_el + labour + maintenance + waste;   CAPEX = f_Lang·Σ C_ref,j·(S/S_ref)ⁿ' },
    { kind: 'formula', label: 'Unit cost, minimum selling price and break-even', text: 'c = (OPEX + CRF·CAPEX)/Q,  CRF = i(1+i)ᴺ/((1+i)ᴺ − 1);   Q_BE = (fixed + CRF·CAPEX)/(p − c_var)' });
  return {
    summary: `Producing ${F(cap, 4)} kg/yr of ${chem.prettyFormula(iv.pr[0])} needs ${F(iv.mPrec, 3)} kg of ${chem.prettyFormula(rt.re[0].buy || iv.re[0])} (${F(pPrec, 3)} USD/kg${indexed ? ', indexed to the metal price' : ''}), ${F(iv.mMake / iv.sv.rho * 1e3, 3)} L of solvent make-up and ${F(iv.eTot, 3)} kWh per kg of product. Operating cost is ${money(opex)} USD/yr and the fixed capital investment ${money(capex)} USD (Lang factor ${LF}); the levelised unit cost, which equals the minimum selling price before tax, is ${F(unit, 4)} USD/kg (${F(unit / inp.fx, 4)} EUR/kg), dominated by ${top[0].toLowerCase()} (${pct(top[1] / unit)}). At ${inp.price} USD/kg the margin is ${pct(margin)} and break-even is reached at ${F(be, 4)} kg/yr.`,
    kpis: [
      kp('unit', 'Levelised unit cost', unit, 'USD/kg', `Annual operating cost plus the capital recovery charge (CRF = ${F(CRF, 3)}) divided by annual output; equal to ${F(unit / inp.fx, 4)} EUR/kg. It is the pre-tax minimum selling price at a ${inp.rate} % discount rate.`, unit < inp.price),
      kp('opex', 'Operating expenditure', opex, 'USD/yr', `Materials, solvent, energy, waste treatment, labour and maintenance at full capacity (${F(opex / cap, 4)} USD/kg), excluding capital charges.`, null),
      kp('capex', 'Fixed capital investment', capex, 'USD', `Purchased equipment cost of ${money(purch)} USD, scaled from the reference capacity with exponent ${inp.scaleExp}, multiplied by the Lang factor ${LF}; a factored estimate with an accuracy of about ±30 %.`, null),
      kp('m_prec', 'Precursor demand', iv.mPrec, 'kg/kg', `Mass of ${chem.prettyFormula(rt.re[0].buy || iv.re[0])} per kilogram of product from the balanced equation at ${pct(iv.y)} yield; the metal content of the precursor is ${pct(wMetal)}.`, null),
      kp('solv', 'Solvent make-up', iv.mMake / iv.sv.rho * 1e3, 'L/kg', `Fresh solvent per kilogram of product after ${inp.solvRec} % recovery of the ${F(iv.Vtot * 1e3, 3)} L used for reaction and washing.`, null),
      kp('energy', 'Process energy', iv.eTot, 'kWh/kg', `Heating of the reaction liquid (${F(iv.eHeat, 3)} kWh), solvent recovery (${F(iv.eRecov, 3)} kWh) and other electricity (${inp.eExtra} kWh) per kilogram of product.`, iv.eTot < 50),
      kp('atom', 'Atom economy', iv.atomEco, '–', 'Mass of product divided by the mass of all reactants in the balanced equation: the stoichiometric ceiling of material efficiency, independent of yield.', iv.atomEco > 0.5),
      kp('efactor', 'E-factor (excluding water)', iv.waste, 'kg/kg', 'Mass of by-product salts, unreacted reagents and unrecovered organic solvent per kilogram of product; aqueous effluent is counted separately.', iv.waste < 5),
      kp('margin', 'Margin at the selling price', margin, '–', `(Price − levelised unit cost)/price at ${inp.price} USD/kg; negative values mean the price does not recover operating cost plus capital.`, margin > 0.15),
      kp('be', 'Break-even volume', be, 'kg/yr', `Annual volume at which the contribution margin (${F(contrib, 3)} USD/kg) covers fixed costs and the capital charge; ${pct(be / cap)} of capacity.`, be > 0 && be < cap),
      kp('m2', 'Material cost per square metre', perM2, 'USD/m²', `Unit cost multiplied by the areal loading of ${inp.loading} g/m².`, null),
      kp('dev', 'Cost per device', perDev, 'USD', `Material cost of a ${inp.devArea} cm² device plus ${inp.integ} USD integration cost.`, null),
    ],
    series: [
      { key: 'brk', title: 'Unit-cost breakdown', xlabel: 'cost item (1 precursor, 2 reagent, 3 solvent, 4 energy, 5 waste, 6 labour, 7 maintenance, 8 capital)', ylabel: 'USD/kg', lines: [{ name: 'cost contribution', x: cats.map((_, i) => i + 1), y: cats.map(c => c[1]) }], explain: 'Contribution of each cost item to the levelised unit cost. The tallest item is where a price negotiation, a yield gain or a design change is most effective.' },
      { key: 'scale', title: 'Unit cost versus plant capacity', xlabel: 'capacity (kg/yr)', ylabel: 'USD/kg', logx: true, logy: true, lines: [{ name: 'levelised unit cost', x: caps, y: caps.map(unitAt) }, { name: 'selling price', x: [caps[0], caps.at(-1)], y: [inp.price, inp.price] }, { name: 'variable cost', x: [caps[0], caps.at(-1)], y: [Math.max(varKg, 1e-9), Math.max(varKg, 1e-9)] }], explain: 'Economies of scale: fixed and capital costs per kilogram fall with capacity while variable costs stay constant. The crossing with the selling price is the smallest viable plant.' },
      { key: 'metal', title: 'Unit cost versus precursor price', xlabel: 'precursor (metal) price relative to today', ylabel: 'USD/kg', lines: [{ name: 'levelised unit cost', x: ms, y: ms.map(unitMetal) }, { name: 'selling price', x: [0.5, 2], y: [inp.price, inp.price] }], explain: 'Exposure of the unit cost to the precursor price, which follows the metal market when indexed. A steep line signals that supply contracts or hedging matter more than process improvements.' },
    ],
    tables: [
      { key: 'mass', title: 'Stoichiometric mass balance per kg of product', columns: ['Species', 'Role', 'Coefficient', 'Mass (kg per kg product)'], rows: iv.rows.map(r => [r[0], r[1], r[2], F(r[3], 4)]), explain: 'Reactant demand and by-product formation per kilogram of saleable product from the balanced equation, including yield and reagent excess. Purchased hydrate forms are counted at their full molar mass.' },
      { key: 'cost', title: 'Cost breakdown', columns: ['Item', 'USD/kg', 'USD/yr', 'Share (%)'], rows: cats.map(([nm, v]) => [nm, F(v, 4), F(v * cap, 4), F(100 * v / unit, 3)]), explain: 'Annual and specific cost by category with its share of the levelised unit cost.' },
      { key: 'eq', title: 'Equipment list and capital estimate', columns: ['Item', 'Cost at reference capacity (USD)', 'Scaled cost (USD)'], rows: [...eq.map(e => [e[0], F(e[1], 4), F(e[2], 4)]), ['Purchased equipment total', F(N.sum(eq.map(e => e[1])), 4), F(purch, 4)], [`Fixed capital (× Lang ${LF})`, '–', F(capex, 4)]], explain: 'Purchased equipment split by typical shares and scaled with the cost–capacity exponent, then multiplied by the Lang factor to give installed fixed capital.' },
    ],
    equations: eqs,
    checks: [
      { name: 'Reaction mass balance', pass: iv.bal ? chem.auditEquation(iv.re, iv.pr, iv.bal.coefficients).massBalanced : true, detail: iv.bal ? `coefficients ${iv.bal.coefficients.join(' : ')}` : 'physical route: feedstock utilisation only' },
      { name: 'Cost items sum to the unit cost', pass: Math.abs(N.sum(cats.map(c => c[1])) - unit) < 1e-9 * unit, detail: `Σ = ${F(N.sum(cats.map(c => c[1])), 6)} USD/kg` },
      { name: 'Selling price covers variable cost', pass: contrib > 0, detail: `contribution margin ${F(contrib, 4)} USD/kg` },
      { name: 'Break-even below capacity', pass: be > 0 && be <= cap, detail: `${F(be, 4)} vs ${F(cap, 4)} kg/yr` },
      { name: 'Capacity within the range of the scaling law', pass: cap / inp.capRef <= 10 && cap / inp.capRef >= 0.1, detail: `capacity ratio ${F(cap / inp.capRef, 3)}` },
    ],
    warnings: w,
    recommendations: [
      { priority: 'high', text: `Attack the largest cost item, ${top[0].toLowerCase()} (${pct(top[1] / unit)} of unit cost): ${top[0] === 'Precursor' ? `each percentage point of yield is worth ${F(iv.mPrec * pPrec * 0.01, 3)} USD/kg, and a cheaper zinc, copper or silver source (oxide, recycled metal, recovered mother liquor) lowers the premium over metal value` : top[0] === 'Capital charge' || top[0] === 'Labour and overhead' ? `it falls with throughput, so raise the precursor concentration or shorten the cycle before adding equipment` : top[0] === 'Energy' ? 'recover heat from the discharged liquid and raise the solids concentration' : 'review the recipe and supply contract for this item'}.`, basis: 'Unit-cost breakdown.' },
      { priority: margin < 0.15 ? 'high' : 'medium', text: margin < 0.15 ? `The margin is ${pct(margin)}: the plant needs ${F(be, 4)} kg/yr to break even. Either a price above ${F(unit * 1.2, 3)} USD/kg or a capacity where unit cost falls below the price is required before investing.` : `The margin of ${pct(margin)} tolerates a precursor price rise of ${pct((inp.price - unit) / Math.max(iv.mPrec * pPrec, 1e-9))} before the price no longer covers the full cost.`, basis: 'Break-even and price sensitivity.' },
      { priority: 'medium', text: `Raise solvent recovery from ${inp.solvRec} % towards 90 % and the precursor concentration above ${F(inp.cPrec, 3)} mol/m³: together they cut the ${F(iv.Vtot * 1e3, 3)} L of solvent handled per kilogram, the heating energy of ${F(iv.eHeat, 3)} kWh/kg and the effluent volume.`, basis: 'Solvent and energy inventory per kilogram of product.', sustainability: true },
      { priority: 'low', text: 'Replace the factored capital estimate by vendor quotations for the main items before a decision; factored estimates carry about ±30 % uncertainty, which the discounted-cash-flow mode can propagate.', basis: 'Accuracy class of Lang-factor estimates.' },
    ],
    publish: { 'econ.capex': capex, 'econ.opex': opex, 'fab.cost_per_kg': unit, 'econ.msp': unit, 'fab.energy_kWh_kg': iv.eTot, 'econ.price': inp.price, 'material.formula': iv.pr[0], 'material.name': rt.name },
    numerics: { method: 'Algebraic mass and energy balance from the balanced reaction (chem kernel); factored capital estimate (six-tenths rule and Lang factor); annuity capital recovery', dof: cats.length, steps: 1, iterations: 1, residual: Math.abs(N.sum(cats.map(c => c[1])) - unit), notes: 'Steady full-capacity operation; prices constant; no tax, working capital or ramp-up (see the cash-flow mode).' },
  };
}

// ---------- mode: discounted cash flow ----------
function cashFlows(P) {
  const cf = [-P.capex]; let wcPrev = 0, dep = 0;
  for (let t = 1; t <= P.life; t++) { const u = P.ramp > 0 ? Math.min(1, t / P.ramp) : 1, rev = P.price * P.vol * u, op = P.opex * ((1 - P.vf) + P.vf * u * P.vol / P.vol0), d = t <= P.depY ? P.capex / P.depY : 0, ebt = rev - op - d, taxp = ebt > 0 ? ebt * P.tax : 0, wcN = P.wc * rev; let c = rev - op - taxp - (wcN - wcPrev); wcPrev = wcN; if (t === P.life) c += wcN; cf.push(c); dep += d; }
  return cf;
}
const npvOf = (cf, r) => cf.reduce((s, c, t) => s + c / Math.pow(1 + r, t), 0);
function runDcf(inp) {
  const w = [], life = Math.round(N.clamp(inp.life, 2, 40)), r = inp.rate / 100, P0 = { capex: inp.capex, price: inp.price, vol: inp.cap, vol0: inp.cap, opex: inp.opex, vf: inp.varFrac / 100, ramp: Math.round(inp.ramp), life, tax: inp.tax / 100, depY: Math.round(inp.depYears), wc: inp.wc / 100 };
  const cf = cashFlows(P0), npv = npvOf(cf, r); let irr = null; const lo = npvOf(cf, -0.9), hi = npvOf(cf, 10);
  if (lo * hi < 0) { try { irr = N.brent(x => npvOf(cf, x), -0.9, 10, 1e-12); } catch { irr = null; } }
  if (irr === null) w.push('The cash flows have no sign change in net present value between −90 % and +1000 % discount rate: the internal rate of return is not defined (the investment is never recovered, or never negative).');
  let cum = 0, dpb = null, spb = null, cumU = 0; const cumD = [], cumS = []; cf.forEach((c, t) => { const d = c / Math.pow(1 + r, t), prev = cum, prevU = cumU; cum += d; cumU += c; cumD.push(cum); cumS.push(cumU); if (dpb === null && t > 0 && cum >= 0) dpb = t - 1 + (d !== 0 ? -prev / d : 0); if (spb === null && t > 0 && cumU >= 0) spb = t - 1 + (c !== 0 ? -prevU / c : 0); });
  const pi = 1 + npv / inp.capex, annuity = npv * crf(r, life), resid = irr !== null ? npvOf(cf, irr) : 0;
  if (P0.depY > life) w.push(`Depreciation period (${P0.depY} yr) exceeds the project life: ${money(inp.capex * (1 - life / P0.depY))} USD remain undepreciated and their tax shield is lost in this model.`);
  // Monte Carlo (mean-preserving lognormal multipliers)
  const rnd = N.rng(Math.round(inp.seed)), nM = Math.round(N.clamp(inp.nMc, 200, 30000)), sg = cv => Math.sqrt(Math.log(1 + (cv / 100) ** 2)), mult = s => s > 0 ? Math.exp(s * rnd.normal() - s * s / 2) : 1, S = [sg(inp.cvPrice), sg(inp.cvVol), sg(inp.cvOpex), sg(inp.cvCapex)], mc = new Array(nM);
  for (let i = 0; i < nM; i++) mc[i] = npvOf(cashFlows({ ...P0, price: P0.price * mult(S[0]), vol: P0.vol * mult(S[1]), opex: P0.opex * mult(S[2]), capex: P0.capex * mult(S[3]) }), r);
  const pNeg = mc.filter(v => v < 0).length / nM, q = [0.05, 0.5, 0.95].map(p => N.quantile(mc, p)), mMean = N.mean(mc), mSd = N.std(mc), h = N.histogram(mc, 30);
  // tornado
  const sw = inp.swing / 100, vars = [['Selling price', 'price'], ['Sales volume', 'vol'], ['Operating cost', 'opex'], ['Capital investment', 'capex'], ['Discount rate', 'r'], ['Tax rate', 'tax']];
  const tor = vars.map(([nm, k]) => { const f = m => k === 'r' ? npvOf(cf, r * m) : npvOf(cashFlows({ ...P0, [k]: k === 'tax' ? Math.min(P0.tax * m, 0.95) : P0[k] * m }), r), a = f(1 - sw), b = f(1 + sw); return { nm, lo: a, hi: b, range: Math.abs(b - a) }; }).sort((a, b) => b.range - a.range);
  const rs = N.linspace(0, Math.max(0.4, (irr || 0) * 1.6), 50), yrs = cf.map((_, t) => t), depSum = Math.min(life, P0.depY) * inp.capex / P0.depY;
  // price for NPV = 0
  let pBE = inp.price; try { pBE = N.brent(p => npvOf(cashFlows({ ...P0, price: p }), r), 1e-6, inp.price * 1000); } catch { pBE = 0; }
  return {
    summary: `Over ${life} years at a ${inp.rate} % discount rate the project has a net present value of ${money(npv)} USD (${money(npv / inp.fx)} EUR), ${irr !== null ? `an internal rate of return of ${pct(irr)}` : 'no defined internal rate of return'}, a discounted payback of ${dpb !== null ? F(dpb, 3) + ' years' : 'more than the project life'} and a profitability index of ${F(pi, 3)}. With the stated uncertainties the Monte Carlo analysis gives a ${pct(pNeg)} probability of a negative NPV (5th–95th percentile ${money(q[0])} to ${money(q[2])} USD). The most influential input is ${tor[0].nm.toLowerCase()}: ±${inp.swing} % moves the NPV between ${money(Math.min(tor[0].lo, tor[0].hi))} and ${money(Math.max(tor[0].lo, tor[0].hi))} USD.`,
    kpis: [
      kp('npv', 'Net present value', npv, 'USD', `Sum of all after-tax cash flows discounted at ${inp.rate} %, including ramp-up, straight-line depreciation over ${P0.depY} years and working capital; positive values create value above the required return.`, npv > 0),
      kp('irr', 'Internal rate of return', irr !== null ? 100 * irr : 'not defined', irr !== null ? '%' : '', irr !== null ? `Discount rate at which the net present value is zero (Brent root search, residual ${F(resid, 2)} USD); to be compared with the required return of ${inp.rate} %.` : 'No discount rate between −90 % and +1000 % makes the net present value zero, so the internal rate of return does not exist for these cash flows.', irr !== null ? irr > r : false),
      kp('dpb', 'Discounted payback', dpb !== null ? dpb : 'not reached', dpb !== null ? 'yr' : '', dpb !== null ? `Time at which the cumulative discounted cash flow returns to zero (linear interpolation within the year); the undiscounted payback is ${spb !== null ? F(spb, 3) + ' years' : 'not reached'}.` : 'The cumulative discounted cash flow stays negative for the whole project life: the investment is not recovered at this discount rate.', dpb !== null && dpb < life * 0.6),
      kp('pi', 'Profitability index', pi, '–', 'Present value of the cash inflows divided by the capital investment (1 + NPV/CAPEX); values above 1 are acceptable and the index ranks projects competing for limited capital.', pi > 1),
      kp('annuity', 'Equivalent annual value', annuity, 'USD/yr', 'Net present value expressed as a constant annual amount over the project life (NPV × capital recovery factor).', annuity > 0),
      kp('p_be', 'Price for zero NPV', pBE, 'USD/kg', `Selling price at which the after-tax net present value is exactly zero: the minimum selling price including tax, ramp-up and working capital (present price ${inp.price} USD/kg).`, pBE > 0 && pBE < inp.price),
      kp('p_neg', 'Probability of negative NPV', pNeg, '–', `Share of ${nM} Monte Carlo scenarios with NPV below zero for lognormal uncertainty of price, volume, operating and capital cost.`, pNeg < 0.2),
      kp('mc_p5', 'NPV, 5th percentile', q[0], 'USD', 'Net present value exceeded in 95 % of the simulated scenarios: a downside figure for risk-averse decisions.', q[0] > 0),
      kp('mc_p50', 'NPV, median', q[1], 'USD', `Median of the simulated NPV distribution (mean ${money(mMean)} USD, standard deviation ${money(mSd)} USD).`, q[1] > 0),
      kp('mc_p95', 'NPV, 95th percentile', q[2], 'USD', 'Net present value exceeded in only 5 % of the scenarios: the upside of the distribution.', null),
      kp('tor', `NPV swing from ${tor[0].nm.toLowerCase()}`, tor[0].range, 'USD', `Change of the NPV when ${tor[0].nm.toLowerCase()} alone varies by ±${inp.swing} %; the largest bar of the tornado ranking.`, null),
    ],
    series: [
      { key: 'cum', title: 'Cumulative cash flow', xlabel: 'year', ylabel: 'USD', lines: [{ name: 'discounted', x: yrs, y: cumD }, { name: 'undiscounted', x: yrs, y: cumS }, { name: 'zero', x: [0, life], y: [0, 0] }], explain: 'Cumulative after-tax cash flow with and without discounting. The crossing of zero is the payback time and the end value of the discounted curve is the net present value.' },
      { key: 'npvr', title: 'Net present value versus discount rate', xlabel: 'discount rate (–)', ylabel: 'NPV (USD)', lines: [{ name: 'NPV', x: rs, y: rs.map(x => npvOf(cf, x)) }, { name: 'zero', x: [rs[0], rs.at(-1)], y: [0, 0] }, ...(irr !== null ? [{ name: 'internal rate of return', x: [irr], y: [0] }] : [])], explain: 'NPV profile: the rate at which the curve crosses zero is the internal rate of return. A steep profile means the result depends strongly on the cost of capital.' },
      { key: 'hist', title: 'Monte Carlo distribution of NPV', xlabel: 'NPV (USD)', ylabel: 'scenarios', lines: [{ name: 'NPV histogram', x: h.x, y: h.y }, { name: 'zero NPV', x: [0, 0], y: [0, N.maxOf(h.y)] }], explain: 'Distribution of the net present value under joint uncertainty of price, volume and costs. The area left of zero is the probability of destroying value.' },
      { key: 'tor', title: 'Tornado sensitivity', xlabel: 'input rank (1 = most influential)', ylabel: 'NPV (USD)', lines: [{ name: `input −${inp.swing} %`, x: tor.map((_, i) => i + 1), y: tor.map(t => t.lo) }, { name: `input +${inp.swing} %`, x: tor.map((_, i) => i + 1), y: tor.map(t => t.hi) }, { name: 'base case', x: [1, tor.length], y: [npv, npv] }], explain: 'NPV when each input is moved down and up by the swing with all others at base values, ranked by spread (names in the tornado table). Wide spreads identify the inputs worth securing by contract or further study.' },
    ],
    tables: [
      { key: 'cf', title: 'Cash-flow projection', columns: ['Year', 'Utilisation', 'Revenue', 'Operating cost', 'After-tax cash flow', 'Discounted', 'Cumulative discounted'], rows: cf.map((c, t) => { const u = t === 0 ? 0 : (P0.ramp > 0 ? Math.min(1, t / P0.ramp) : 1); return [t, +u.toFixed(2), F(P0.price * P0.vol * u, 4), F(t === 0 ? 0 : P0.opex * ((1 - P0.vf) + P0.vf * u), 4), F(c, 4), F(c / Math.pow(1 + r, t), 4), F(cumD[t], 4)]; }), explain: 'Year-by-year revenue, operating cost and after-tax cash flow with discounting. Year 0 carries the investment; the last year includes recovery of working capital.' },
      { key: 'tor', title: 'Tornado sensitivity table', columns: ['Rank', 'Input', `NPV at −${inp.swing} %`, `NPV at +${inp.swing} %`, 'Range'], rows: tor.map((t, i) => [i + 1, t.nm, F(t.lo, 4), F(t.hi, 4), F(t.range, 4)]), explain: 'One-at-a-time sensitivity of the NPV, ranked by range. It shows direction as well as magnitude: revenue drivers raise NPV when increased, cost drivers lower it.' },
    ],
    equations: [
      { kind: 'formula', label: 'Net present value', text: 'NPV = Σₜ₌₀ᴺ CFₜ/(1 + i)ᵗ,   CFₜ = (Rₜ − Oₜ − Dₜ)·(1 − τ) + Dₜ − ΔWCₜ,   CF₀ = −CAPEX' },
      { kind: 'formula', label: 'Internal rate of return and payback', text: 'NPV(IRR) = 0;   discounted payback: Σₜ CFₜ/(1 + i)ᵗ = 0;   PI = 1 + NPV/CAPEX' },
      { kind: 'formula', label: 'Monte Carlo sampling', text: 'x = x̄·exp(σ·z − σ²/2),  σ² = ln(1 + CV²),  z ~ N(0,1);   P(NPV < 0) = #{NPV < 0}/n' },
    ],
    checks: [
      { name: 'IRR reproduces NPV = 0', pass: irr === null || Math.abs(resid) < 1e-6 * Math.max(inp.capex, 1), detail: irr !== null ? `NPV(IRR) = ${F(resid, 2)} USD` : 'IRR not defined' },
      { name: 'Depreciation recovers the investment', pass: Math.abs(depSum - inp.capex) < 1e-6 * inp.capex, detail: `Σ depreciation ${F(depSum, 5)} vs CAPEX ${F(inp.capex, 5)} USD` },
      { name: 'Monte Carlo mean consistent with base case', pass: Math.abs(mMean - npv) < 4 * mSd / Math.sqrt(nM) + 0.1 * Math.abs(npv) + 1e-6 * inp.capex, detail: `mean ${F(mMean, 4)} vs base ${F(npv, 4)} USD (tax asymmetry causes a small bias)` },
      { name: 'Positive net present value', pass: npv > 0, detail: `${F(npv, 4)} USD` },
      { name: 'Probability of loss below 20 %', pass: pNeg < 0.2, detail: pct(pNeg) },
    ],
    warnings: w,
    recommendations: [
      { priority: npv > 0 ? 'medium' : 'high', text: npv > 0 ? `The base case creates value (NPV ${money(npv)} USD${irr !== null ? `, IRR ${pct(irr)}` : ''}); secure ${tor[0].nm.toLowerCase()} first, because a ${inp.swing} % adverse move changes the NPV by ${money(Math.abs(npv - Math.min(tor[0].lo, tor[0].hi)))} USD.` : `The base case destroys value: a selling price of at least ${F(pBE, 4)} USD/kg, or the equivalent reduction in costs, is needed for a zero NPV.`, basis: 'Net present value and tornado ranking.' },
      { priority: pNeg > 0.2 ? 'high' : 'low', text: pNeg > 0.2 ? `A ${pct(pNeg)} chance of loss is high for a single investment: stage the capital (pilot line first), seek an off-take agreement to narrow price and volume uncertainty, or share equipment with existing production.` : `The probability of a negative NPV is ${pct(pNeg)}; the downside (5th percentile) is ${money(q[0])} USD.`, basis: 'Monte Carlo distribution of NPV.' },
      { priority: 'medium', text: `Shorten the ramp-up: each year of delayed full output costs roughly ${money(Math.abs(npv - npvOf(cashFlows({ ...P0, ramp: P0.ramp + 1 }), r)))} USD of NPV, usually more than modest savings in capital.`, basis: 'Sensitivity of discounted cash flow to early-year revenue.' },
      { priority: 'medium', text: 'Include a carbon price and the avoided cost of solvent, effluent and waste when comparing process options: measures that cut emissions and waste (solvent recovery, heat integration, renewable power contracts) often also stabilise the operating cost against energy-price swings.', basis: 'Operating-cost sensitivity in the tornado table.', sustainability: true },
    ],
    publish: { 'econ.npv': npv, ...(irr !== null ? { 'econ.irr': irr } : {}), ...(dpb !== null ? { 'econ.payback_y': dpb } : {}), 'econ.p_loss': pNeg, 'econ.capex': inp.capex, 'econ.opex': inp.opex },
    numerics: { method: 'Annual discounted cash-flow model; Brent root search for IRR and zero-NPV price; seeded Monte Carlo with lognormal multipliers; one-at-a-time tornado', dof: life + 1, steps: nM, iterations: nM, residual: Math.abs(resid), notes: 'End-of-year cash flows, constant prices (real terms), no loss carry-forward, uncertain inputs independent and constant over the project.' },
  };
}

// ---------- mode: life-cycle screening and circularity ----------
function footprint(inp, over = {}) {
  const q = { ...inp, ...over }, iv = inventory(inp, over), g = q.grid / 1e3, tkm = (iv.mPrec + iv.mReag + iv.mMake * (iv.organic ? 1 : 0)) / 1e3 * q.transKm;
  const st = [['Precursor production', iv.mPrec * q.efPrec], ['Second reagent', iv.mReag * q.efReag], ['Solvent make-up', iv.mMake * q.efSolv], ['Electricity and heat', iv.eTot * g], ['Waste treatment', iv.waste * q.efWaste], ['Transport of inputs', tkm * q.efTrans]];
  const ced = iv.mPrec * q.cedPrec + iv.mReag * q.cedReag + iv.eTot * q.pef + (iv.organic ? iv.mMake * 30 : 0), water = (iv.organic ? 0 : iv.mMake) + q.coolWater;
  return { iv, st, gwp: N.sum(st.map(s => s[1])), ced, water };
}
function runLca(inp) {
  const w = [], base = footprint(inp), iv = base.iv, fuK = inp.fu === 'm2' ? inp.loading / 1e3 : inp.fu === 'device' ? inp.loading / 1e3 * inp.devArea / 1e4 : 1, fuName = { kg: 'kg', m2: 'm²', device: 'device' }[inp.fu] || 'kg';
  const sc = [['Baseline', {}], ['No solvent recovery', { solvRec: 0 }], ['Solvent recovery 90 %', { solvRec: 90 }], ['Renewable electricity', { grid: inp.gridRen }], ['Heat recovery 60 %', { heatRec: 60 }], ['Yield 98 %', { yieldFrac: Math.max(inp.yieldFrac, 0.98) }], ['Combined measures', { solvRec: 90, grid: inp.gridRen, heatRec: 60, yieldFrac: Math.max(inp.yieldFrac, 0.98) }]].map(([nm, o]) => { const f = footprint(inp, o); return { nm, gwp: f.gwp, ced: f.ced, water: f.water }; });
  const top = base.st.slice().sort((a, b) => b[1] - a[1])[0], pay = inp.eSave > 0 && inp.grid > 0 ? base.gwp / (inp.eSave * inp.grid / 1e3) : 0;
  // material circularity indicator (Ellen MacArthur Foundation formulation), per unit mass
  const FR = inp.recIn / 100, CR = inp.recOut / 100, EC = inp.recEff / 100, X = Math.max(inp.lifeRatio, 1e-6), lfiOf = (fr, cr) => { const V = 1 - fr, WC = (1 - EC) * cr, WF = (1 - EC) * fr / EC, W = 1 - cr + (WF + WC) / 2; return N.clamp((V + W) / (2 + (WF - WC) / 2), 0, 1); }, mciOf = (fr, cr) => Math.max(0, 1 - lfiOf(fr, cr) * 0.9 / X), LFI = lfiOf(FR, CR), mci = mciOf(FR, CR), mciIn = mciOf(Math.min(1, FR + 0.3), CR), mciOut = mciOf(FR, Math.min(1, CR + 0.3));
  // critical raw materials
  const els = new Set(); const addF = f => { try { Object.keys(chem.parseFormula(f).atoms).forEach(e => els.add(e)); } catch (e) { w.push(`Formula "${f}" could not be read for the critical-material screening (${e.message}).`); } };
  [iv.pr[0], iv.rt.re[0].buy || iv.re[0]].forEach(addF); if (String(inp.elements || '').trim()) String(inp.elements).split(/[;,\s]+/).filter(Boolean).forEach(addF);
  const flagged = [...els].filter(e => CRM[e]), eqs = routeEquations(iv);
  if (flagged.length) w.push(`Elements on critical-raw-material screening lists: ${flagged.map(e => `${e} (${CRM[e]})`).join(', ')}. Verify against the current official list and assess supply risk and substitution.`);
  w.push('Screening-level cradle-to-gate estimate with user-adjustable emission and energy factors: it identifies hot spots and compares scenarios but is not a life-cycle assessment conforming to ISO 14040/14044 (no critical review, limited impact categories, capital goods excluded).');
  if (top[1] / base.gwp > 0.8) w.push(`${top[0]} carries ${pct(top[1] / base.gwp)} of the footprint, so the result is only as reliable as that single factor: obtain supplier-specific data for it.`);
  eqs.push({ kind: 'formula', label: 'Inventory and impact', text: 'GWP = Σ mᵢ·EFᵢ + E·EF_grid + m_waste·EF_waste + (t·km)·EF_transport;   CED = Σ mᵢ·CEDᵢ + E·PEF' },
    { kind: 'formula', label: 'Material circularity indicator', text: 'V = M(1 − F_R);  W = M(1 − C_R) + (W_F + W_C)/2;  LFI = (V + W)/(2M + (W_F − W_C)/2);  MCI = max(0, 1 − LFI·0.9/X)', note: 'W_C = M(1 − E_C)C_R, W_F = M(1 − E_F)F_R/E_F, X = lifetime or use intensity relative to the industry average; reuse and composting flows are taken as zero.' },
    { kind: 'formula', label: 'Carbon payback', text: 't_payback = GWP_production/(E_saved·EF_grid)' });
  const recs = N.linspace(0, 95, 20), grids = N.linspace(0, Math.max(inp.grid * 1.5, 100), 20);
  return {
    summary: `Cradle-to-gate screening for ${chem.prettyFormula(iv.pr[0])}: ${F(base.gwp * fuK, 3)} kg CO₂e, ${F(base.ced * fuK, 3)} MJ primary energy and ${F(base.water * fuK, 3)} L of water per ${fuName}. The largest contribution is ${top[0].toLowerCase()} (${pct(top[1] / base.gwp)}). Solvent recovery of 90 % changes the footprint to ${F(sc[2].gwp * fuK, 3)}, renewable electricity to ${F(sc[3].gwp * fuK, 3)} and all measures combined to ${F(sc[6].gwp * fuK, 3)} kg CO₂e per ${fuName} (${pct(1 - sc[6].gwp / base.gwp)} lower). The material circularity indicator is ${F(mci, 3)}${flagged.length ? ` and ${flagged.length} element(s) appear on critical-raw-material screening lists` : ''}. Factors are user-adjustable screening values.`,
    kpis: [
      kp('gwp', `Global-warming potential per ${fuName}`, base.gwp * fuK, 'kg CO₂e', `Greenhouse-gas emissions from precursor and reagent production, solvent, energy, waste treatment and transport; ${F(base.gwp, 3)} kg CO₂e per kg of product.`, null),
      kp('ced', `Cumulative energy demand per ${fuName}`, base.ced * fuK, 'MJ', `Primary energy embodied in the inputs plus process electricity converted with a primary energy factor of ${inp.pef} MJ/kWh.`, null),
      kp('water', `Water use per ${fuName}`, base.water * fuK, 'L', `Process-water make-up after ${inp.solvRec} % recovery plus cooling and utility water.`, null),
      kp('top', `Share of ${top[0].toLowerCase()}`, top[1] / base.gwp, '–', 'Contribution of the largest life-cycle stage to the global-warming potential: the hot spot where reduction measures and better data matter most.', top[1] / base.gwp < 0.6),
      kp('gwp_rec', 'GWP with 90 % solvent recovery', sc[2].gwp * fuK, 'kg CO₂e', `Footprint if 90 % of the solvent is recovered (including the energy of recovery); without any recovery it would be ${F(sc[1].gwp * fuK, 3)} kg CO₂e.`, sc[2].gwp <= base.gwp),
      kp('gwp_ren', 'GWP with renewable electricity', sc[3].gwp * fuK, 'kg CO₂e', `Footprint with an electricity carbon intensity of ${inp.gridRen} g/kWh instead of ${inp.grid} g/kWh.`, sc[3].gwp <= base.gwp),
      kp('gwp_best', 'GWP with all measures combined', sc[6].gwp * fuK, 'kg CO₂e', 'Footprint with 90 % solvent recovery, 60 % heat recovery, renewable electricity and at least 98 % yield together.', sc[6].gwp <= base.gwp),
      kp('payback', 'Carbon payback time', pay, 'yr', inp.eSave > 0 && inp.grid > 0 ? `Years of use after which the electricity saved by the material in its application (${inp.eSave} kWh/yr per kg at ${inp.grid} g CO₂/kWh) offsets its production footprint.` : 'No energy saving in use or zero grid intensity entered: no carbon payback applies and zero is shown.', inp.eSave > 0 ? pay < 3 : null),
      kp('mci', 'Material circularity indicator', mci, '–', `Circularity on a 0–1 scale from recycled content (${inp.recIn} %), collection (${inp.recOut} %), recycling efficiency (${inp.recEff} %) and relative lifetime (${inp.lifeRatio}); 0.1 corresponds to a fully linear product of average life.`, mci > 0.5),
      kp('lfi', 'Linear flow index', LFI, '–', 'Share of the material flowing linearly from virgin source to unrecoverable waste (1 = fully linear, 0 = fully restorative).', LFI < 0.5),
      kp('crm', 'Elements on critical-material screening lists', flagged.length, `of ${els.size}`, flagged.length ? `Flagged: ${flagged.join(', ')}. Screening against elements that have appeared on EU critical-raw-material lists; the current official list must be checked.` : 'None of the elements of product, precursor or the additional formula is on the built-in screening list.', flagged.length === 0),
      kp('efactor', 'E-factor (excluding water)', iv.waste, 'kg/kg', 'Mass of by-products, unreacted reagents and unrecovered organic solvent per kilogram of product.', iv.waste < 5),
    ],
    series: [
      { key: 'contrib', title: 'Contribution by life-cycle stage', xlabel: 'stage (1 precursor, 2 reagent, 3 solvent, 4 energy, 5 waste, 6 transport)', ylabel: `kg CO₂e per ${fuName}`, lines: [{ name: 'global-warming potential', x: base.st.map((_, i) => i + 1), y: base.st.map(s => s[1] * fuK) }], explain: 'Greenhouse-gas emissions of each cradle-to-gate stage. Dilute wet-chemical routes are usually dominated by energy and solvent, noble-metal routes by the precursor.' },
      { key: 'scen', title: 'Scenario comparison', xlabel: 'scenario (1 baseline … 7 combined)', ylabel: `kg CO₂e per ${fuName}`, lines: [{ name: 'global-warming potential', x: sc.map((_, i) => i + 1), y: sc.map(s => s.gwp * fuK) }], explain: 'Footprint of the baseline and of each improvement scenario in the order of the scenario table. The drop from scenario 1 to 7 is the improvement potential within the present route.' },
      { key: 'rec', title: 'Footprint versus solvent recovery', xlabel: 'solvent recovery (%)', ylabel: `kg CO₂e per ${fuName}`, lines: [{ name: 'present grid', x: recs, y: recs.map(v => footprint(inp, { solvRec: v }).gwp * fuK) }, { name: 'renewable electricity', x: recs, y: recs.map(v => footprint(inp, { solvRec: v, grid: inp.gridRen }).gwp * fuK) }], explain: 'Effect of solvent recovery on the footprint for two electricity supplies. For organic solvents recovery trades solvent production against distillation energy, so the benefit depends on the grid.' },
      { key: 'grid', title: 'Footprint versus grid carbon intensity', xlabel: 'grid carbon intensity (g CO₂/kWh)', ylabel: `kg CO₂e per ${fuName}`, lines: [{ name: 'global-warming potential', x: grids, y: grids.map(v => footprint(inp, { grid: v }).gwp * fuK) }], explain: 'Sensitivity of the footprint to the electricity mix. The intercept at zero is the part that only material efficiency and precursor choice can reduce.' },
    ],
    tables: [
      { key: 'inv', title: 'Inventory per kg of product', columns: ['Flow', 'Amount', 'Unit'], rows: [['Precursor (' + chem.prettyFormula(iv.rt.re[0].buy || iv.re[0]) + ')', F(iv.mPrec, 4), 'kg'], ['Second reagent', F(iv.mReag, 4), 'kg'], ['Solvent used', F(iv.mSolv, 4), 'kg'], ['Solvent make-up', F(iv.mMake, 4), 'kg'], ['Electricity and heat', F(iv.eTot, 4), 'kWh'], ['By-products and residues', F(iv.waste, 4), 'kg'], ['Aqueous effluent', F(iv.ww * 1e3, 4), 'L']], explain: 'Material and energy flows per kilogram of product from the stoichiometric balance, which the emission factors convert into impacts.' },
      { key: 'scen', title: 'Scenarios', columns: ['#', 'Scenario', `GWP (kg CO₂e/${fuName})`, `CED (MJ/${fuName})`, `Water (L/${fuName})`, 'GWP change (%)'], rows: sc.map((s, i) => [i + 1, s.nm, F(s.gwp * fuK, 4), F(s.ced * fuK, 4), F(s.water * fuK, 4), F(100 * (s.gwp / base.gwp - 1), 3)]), explain: 'Global-warming potential, cumulative energy demand and water use of each scenario relative to the baseline.' },
      { key: 'crm', title: 'Critical-raw-material screening', columns: ['Element', 'On screening list', 'Listed as'], rows: [...els].sort().map(e => [e, CRM[e] ? 'yes' : 'no', CRM[e] || '–']), explain: 'Elements of the product, precursor and any additional formula checked against elements that have appeared on EU critical-raw-material lists. The list changes with each revision, so the official current version governs.' },
    ],
    equations: eqs,
    checks: [
      { name: 'Stage contributions sum to the total', pass: Math.abs(N.sum(base.st.map(s => s[1])) - base.gwp) < 1e-9 * Math.max(base.gwp, 1e-300), detail: `Σ = ${F(base.gwp, 6)} kg CO₂e/kg` },
      { name: 'Circularity indicator within [0, 1]', pass: mci >= 0 && mci <= 1, detail: `MCI = ${F(mci, 4)}, LFI = ${F(LFI, 4)}` },
      { name: 'Combined measures reduce the footprint', pass: sc[6].gwp <= base.gwp * (1 + 1e-9), detail: `${F(sc[6].gwp, 4)} vs ${F(base.gwp, 4)} kg CO₂e/kg` },
      { name: 'No single stage above 80 % (data-quality concentration)', pass: top[1] / base.gwp <= 0.8, detail: `${top[0]}: ${pct(top[1] / base.gwp)}` },
      { name: 'No critical raw materials flagged', pass: flagged.length === 0, detail: flagged.length ? flagged.join(', ') : 'none' },
    ],
    warnings: w,
    recommendations: [
      { priority: 'high', text: `Reduce the hot spot, ${top[0].toLowerCase()} (${pct(top[1] / base.gwp)}): ${top[0].startsWith('Electricity') ? `process energy is ${F(iv.eTot, 3)} kWh/kg, so raise the solids concentration, recover heat and contract low-carbon electricity` : top[0].startsWith('Precursor') ? 'request supplier-specific footprint data, raise the yield and consider recycled or secondary metal feedstock' : top[0].startsWith('Solvent') ? 'increase solvent recovery or switch to water' : 'review this stage for reduction options'}.`, basis: 'Contribution analysis of the global-warming potential.', sustainability: true },
      { priority: 'medium', text: `The combined scenario lowers the footprint by ${pct(1 - sc[6].gwp / base.gwp)}; of the single measures, "${sc.slice(1, -1).sort((a, b) => a.gwp - b.gwp)[0].nm}" is the most effective and should be costed first in the techno-economic mode.`, basis: 'Scenario comparison.', sustainability: true },
      { priority: flagged.length ? 'high' : 'low', text: flagged.length ? `Document sourcing and end-of-life recovery for ${flagged.join(', ')}, evaluate substitution or reduced loading, and include the recovery route in the material passport.` : 'No element of the declared composition is on the screening list; repeat the check when dopants, electrodes or substrates are added.', basis: 'Critical-raw-material screening of the element list.', sustainability: true },
      { priority: 'medium', text: mciOut >= mciIn ? `Circularity (MCI ${F(mci, 3)}) responds most to end-of-life collection: raising it from ${inp.recOut} to ${Math.min(100, inp.recOut + 30)} % lifts the indicator to ${F(mciOut, 3)}, against ${F(mciIn, 3)} for the same increase in recycled content.` : `Circularity (MCI ${F(mci, 3)}) responds most to recycled feedstock: raising recycled content from ${inp.recIn} to ${Math.min(100, inp.recIn + 30)} % lifts the indicator to ${F(mciIn, 3)}, against ${F(mciOut, 3)} for the same increase in collection.`, basis: 'Material circularity indicator and linear flow index.', sustainability: true },
      { priority: 'low', text: 'Before any external environmental claim, commission a life-cycle assessment with qualified inventory data, further impact categories and critical review; this screening is for internal design decisions.', basis: 'Scope limits of the screening estimate.' },
    ],
    publish: { 'econ.co2_per_kg': base.gwp, 'lca.ced_MJ_kg': base.ced, 'lca.water_L_kg': base.water, 'lca.mci': mci, 'fab.energy_kWh_kg': iv.eTot, 'material.formula': iv.pr[0] },
    numerics: { method: 'Process-based inventory from the balanced reaction multiplied by screening emission and energy factors; scenario re-evaluation; closed-form circularity indicator', dof: base.st.length, steps: sc.length, iterations: sc.length + 40, residual: 0, notes: 'Cradle-to-gate system boundary; capital goods, packaging and use phase excluded; single impact category for climate plus energy and water inventories.' },
  };
}

// ---------- mode: learning curve, scale economics, make-or-buy ----------
function runLearning(inp) {
  const w = [], lr = inp.lr / 100, b = Math.log2(Math.max(1 - lr, 1e-9)), c0 = inp.c0, q0 = inp.q0, cQ = Q => c0 * Math.pow(Math.max(Q, 1e-300) / q0, b), cT = cQ(inp.qT);
  const reach = lr > 0 && inp.cTarget < c0, Qs = reach ? q0 * Math.pow(inp.cTarget / c0, 1 / b) : (inp.cTarget >= c0 ? q0 : 0), g = inp.growth / 100, dQ = Math.max(Qs - q0, 0);
  const tReach = !reach ? 0 : g > 1e-9 ? Math.log(1 + g * dQ / inp.qYear) / Math.log(1 + g) : dQ / inp.qYear, doublings = reach ? Math.log2(Qs / q0) : 0;
  const avg = inp.qT > q0 * (1 + 1e-9) ? (Math.abs(1 + b) > 1e-9 ? c0 * q0 / (1 + b) * (Math.pow(inp.qT / q0, 1 + b) - 1) / (inp.qT - q0) : c0 * q0 * Math.log(inp.qT / q0) / (inp.qT - q0)) : c0;
  if (lr === 0) w.push('Learning rate is zero: cost does not fall with cumulative volume and the target cost cannot be reached by experience alone.');
  if (inp.cTarget >= c0) w.push('The target cost is not below the present cost: it is already met.');
  if (reach && tReach > 30) w.push(`At ${inp.growth} %/yr growth the target cost is reached only after ${F(tReach, 3)} years: learning alone will not deliver it within a planning horizon.`);
  if (inp.qT < q0) w.push('Target cumulative production is below the production to date; the curve is evaluated backwards in time.');
  const n = inp.scaleExp, S = inp.sRatio, cs = inp.capShare / 100, capRatio = Math.pow(S, n - 1), cScale = c0 * ((1 - cs) + cs * capRatio);
  const dm = inp.buyPrice - inp.varMake, qBE = dm > 0 ? inp.fixedMake / dm : 0, saving = dm * inp.demandMB - inp.fixedMake, make = saving > 0, costMake = inp.fixedMake / inp.demandMB + inp.varMake;
  if (dm <= 0) w.push('The purchase price is not above the in-house variable cost: buying is cheaper at every volume and no break-even exists (reported as 0).');
  const Qg = N.logspace(Math.log10(Math.min(q0, inp.qT) / 2), Math.log10(Math.max(inp.qT, q0 * 4) * 2), 60), yrs = N.linspace(0, 20, 41), cumAt = t => q0 + (g > 1e-9 ? inp.qYear * (Math.pow(1 + g, t) - 1) / g : inp.qYear * t), Ss = N.logspace(-1, 3, 50), Dm = N.logspace(Math.log10(inp.demandMB / 20), Math.log10(inp.demandMB * 20), 50);
  return {
    summary: `With a learning rate of ${inp.lr} % per doubling (exponent b = ${F(b, 4)}) the unit cost falls from ${F(c0, 4)} USD/kg at ${F(q0, 3)} kg to ${F(cT, 4)} USD/kg at ${F(inp.qT, 3)} kg cumulative production. The target of ${F(inp.cTarget, 4)} USD/kg ${reach ? `needs ${F(Qs, 3)} kg (${F(doublings, 3)} doublings), reached after about ${F(tReach, 3)} years at ${inp.growth} %/yr growth` : inp.cTarget >= c0 ? 'is already met' : 'is not reachable by learning'}. A ${F(S, 3)}-fold larger plant lowers capital cost per kilogram to ${pct(capRatio)} and the unit cost to ${F(cScale, 4)} USD/kg. In-house production is cheaper than buying above ${F(qBE, 4)} kg/yr; at ${F(inp.demandMB, 4)} kg/yr the better choice is to ${make ? 'make' : 'buy'}.`,
    kpis: [
      kp('b', 'Learning exponent b', b, '–', `Exponent of Wright's law, log₂(1 − learning rate); each doubling of cumulative output multiplies cost by ${F(1 - lr, 3)}.`, null),
      kp('c_target_q', 'Unit cost at target volume', cT, 'USD/kg', `Projected cost at ${F(inp.qT, 3)} kg cumulative production, ${pct(1 - cT / c0)} below the present cost. Extrapolation over many doublings is uncertain because learning rates are not constant.`, cT <= inp.cTarget),
      kp('q_star', 'Cumulative volume for target cost', Qs, 'kg', reach ? `Cumulative production at which the experience curve reaches ${F(inp.cTarget, 4)} USD/kg: ${F(doublings, 3)} doublings from today.` : 'The target cost is already met or cannot be reached by learning; the present cumulative volume or zero is shown.', reach ? Qs <= inp.qT : inp.cTarget >= c0),
      kp('t_star', 'Years to reach target cost', tReach, 'yr', reach ? `Calendar time for the cumulative volume to grow from ${F(q0, 3)} to ${F(Qs, 3)} kg with annual output starting at ${F(inp.qYear, 3)} kg/yr and growing ${inp.growth} % per year.` : 'No waiting time applies because the target is already met or unreachable; zero is shown.', reach ? tReach < 10 : inp.cTarget >= c0),
      kp('avg', 'Average cost of the next production', avg, 'USD/kg', `Mean unit cost of everything produced between ${F(q0, 3)} and ${F(inp.qT, 3)} kg (integral of the experience curve): the figure to budget the scale-up programme with.`, null),
      kp('cap_ratio', 'Capital cost per kg at the larger scale', capRatio, '–', `Specific capital cost relative to today for a ${F(S, 3)}-fold capacity with the cost–capacity exponent ${n}: S^(n−1).`, capRatio < 1),
      kp('c_scale', 'Unit cost at the larger scale', cScale, 'USD/kg', `Present unit cost with its capital-related share (${inp.capShare} %) reduced by the scale effect; variable costs per kilogram are unchanged.`, cScale < c0),
      kp('q_be', 'Make-or-buy break-even volume', qBE, 'kg/yr', dm > 0 ? `Annual requirement above which in-house production (fixed ${money(inp.fixedMake)} USD/yr plus ${inp.varMake} USD/kg) is cheaper than buying at ${inp.buyPrice} USD/kg.` : 'Buying is cheaper at every volume because the purchase price does not exceed the in-house variable cost; zero is shown.', dm > 0 && inp.demandMB > qBE),
      kp('saving', 'Annual saving of making over buying', saving, 'USD/yr', `Difference between buying and making ${F(inp.demandMB, 4)} kg/yr; in-house unit cost is ${F(costMake, 4)} USD/kg. Negative values favour buying.`, saving > 0),
      kp('decision', 'Make-or-buy indication', make ? 'make in-house' : 'buy', '', 'Cost-based indication at the stated annual requirement; supply security, quality control, intellectual property and capital availability must be weighed separately.', null),
    ],
    series: [
      { key: 'exp', title: 'Experience curve', xlabel: 'cumulative production (kg)', ylabel: 'unit cost (USD/kg)', logx: true, logy: true, lines: [{ name: `learning rate ${inp.lr} %`, x: Qg, y: Qg.map(cQ) }, { name: 'target cost', x: [Qg[0], Qg.at(-1)], y: [inp.cTarget, inp.cTarget] }, { name: 'today', x: [q0], y: [c0] }], explain: 'Unit cost against cumulative production on logarithmic axes, where Wright\'s law is a straight line. The crossing with the target cost gives the production volume the market must absorb first.' },
      { key: 'time', title: 'Projected unit cost over time', xlabel: 'years from now', ylabel: 'unit cost (USD/kg)', lines: [{ name: 'projected cost', x: yrs, y: yrs.map(t => cQ(cumAt(t))) }, { name: 'target cost', x: [0, 20], y: [inp.cTarget, inp.cTarget] }], explain: 'The experience curve translated into calendar time with the assumed growth of annual output. Fast growth, not time itself, drives the cost down.' },
      { key: 'scale', title: 'Unit cost versus plant scale', xlabel: 'capacity relative to present', ylabel: 'unit cost (USD/kg)', logx: true, lines: [{ name: 'unit cost', x: Ss, y: Ss.map(s => c0 * ((1 - cs) + cs * Math.pow(s, n - 1))) }, { name: 'variable-cost floor', x: [0.1, 1000], y: [c0 * (1 - cs), c0 * (1 - cs)] }], explain: 'Economies of scale from the cost–capacity law: the capital-related part of the cost falls with capacity and the unit cost approaches the variable-cost floor.' },
      { key: 'mb', title: 'Make versus buy', xlabel: 'annual requirement (kg/yr)', ylabel: 'unit cost (USD/kg)', logx: true, logy: true, lines: [{ name: 'make in-house', x: Dm, y: Dm.map(d => inp.fixedMake / d + inp.varMake) }, { name: 'buy', x: [Dm[0], Dm.at(-1)], y: [inp.buyPrice, inp.buyPrice] }], explain: 'In-house unit cost against annual requirement compared with the purchase price. To the right of the crossing the fixed cost is spread thinly enough for in-house production to win.' },
    ],
    tables: [{ key: 'dbl', title: 'Cost by doubling of cumulative production', columns: ['Doublings', 'Cumulative production (kg)', 'Unit cost (USD/kg)', 'Reduction (%)', 'Years from now'], rows: [0, 1, 2, 3, 4, 5, 6, 8, 10].map(k => { const Q = q0 * Math.pow(2, k), dq = Q - q0; return [k, F(Q, 4), F(cQ(Q), 4), F(100 * (1 - cQ(Q) / c0), 3), F(g > 1e-9 ? Math.log(1 + g * dq / inp.qYear) / Math.log(1 + g) : dq / inp.qYear, 3)]; }), explain: 'Projected cost after successive doublings of cumulative production, with the time needed at the assumed growth rate.' }],
    equations: [
      { kind: 'formula', label: "Wright's law", text: 'C(Q) = C₀·(Q/Q₀)ᵇ,  b = log₂(1 − LR);   Q* = Q₀·(C_target/C₀)^(1/b)' },
      { kind: 'formula', label: 'Cumulative production with growth', text: 'Q(t) = Q₀ + q·[(1 + g)ᵗ − 1]/g   ⇒   t* = ln(1 + g·(Q* − Q₀)/q)/ln(1 + g)' },
      { kind: 'formula', label: 'Scale economics', text: 'CAPEX ∝ Sⁿ  ⇒  capital cost per unit ∝ S^(n−1);   c(S) = c₀·[(1 − s_cap) + s_cap·S^(n−1)]' },
      { kind: 'formula', label: 'Make-or-buy break-even', text: 'Q_BE = F/(p_buy − v);   saving = (p_buy − v)·D − F' },
    ],
    checks: [
      { name: 'Target cost reachable by learning', pass: reach || inp.cTarget >= c0, detail: reach ? `${F(doublings, 3)} doublings needed` : inp.cTarget >= c0 ? 'already met' : 'learning rate is zero' },
      { name: 'Target reached within ten years', pass: (reach && tReach <= 10) || inp.cTarget >= c0, detail: reach ? `${F(tReach, 3)} years` : 'see above' },
      { name: 'Experience curve passes through the present point', pass: Math.abs(cQ(q0) - c0) < 1e-9 * c0, detail: `C(Q₀) = ${F(cQ(q0), 6)} USD/kg` },
      { name: 'Scale exponent below 1 (economies of scale exist)', pass: n < 1, detail: `n = ${n}` },
      { name: 'Requirement above make-or-buy break-even', pass: dm > 0 && inp.demandMB > qBE, detail: `${F(inp.demandMB, 4)} vs ${F(qBE, 4)} kg/yr` },
    ],
    warnings: w,
    recommendations: [
      { priority: reach && tReach > 10 ? 'high' : 'medium', text: reach ? `Reaching ${F(inp.cTarget, 4)} USD/kg needs ${F(doublings, 3)} doublings of cumulative output; combine learning with a step change (${F(S, 3)}-fold scale lowers cost to ${F(cScale, 4)} USD/kg at once) rather than waiting ${F(tReach, 3)} years for experience alone.` : 'The target cost is already met at present volume; use the experience curve to set price expectations for larger contracts.', basis: "Wright's law projection and scale economics." },
      { priority: 'medium', text: make ? `At ${F(inp.demandMB, 4)} kg/yr making saves ${money(saving)} USD/yr; confirm that the requirement is secure for the pay-back period, since below ${F(qBE, 4)} kg/yr the fixed cost makes buying cheaper.` : `At ${F(inp.demandMB, 4)} kg/yr buying is cheaper by ${money(-saving)} USD/yr; revisit the decision when the requirement approaches ${F(qBE, 4)} kg/yr or if the supplier cannot guarantee quality and supply.`, basis: 'Make-or-buy break-even analysis.' },
      { priority: 'low', text: `Do not extrapolate the learning rate beyond about five doublings without evidence: mature process steps (precursor supply, energy) set a floor near the variable cost of ${F(c0 * (1 - cs), 4)} USD/kg.`, basis: 'Variable-cost floor of the scale analysis.' },
      { priority: 'medium', text: 'Direct the learning effort at yield, solvent recycling and energy per kilogram: these are the experience gains that lower cost and environmental burden together, whereas savings from cheaper but dirtier inputs do not survive carbon pricing.', basis: 'Sources of cost reduction along the experience curve.', sustainability: true },
    ],
    publish: { 'fab.cost_per_kg': c0, 'econ.cost_at_target': cT, 'econ.learning_b': b, 'econ.make_buy_breakeven': qBE },
    numerics: { method: "Closed-form Wright's-law experience curve, cost–capacity scaling and linear break-even analysis", dof: 1, steps: 60, iterations: 1, residual: Math.abs(cQ(q0) - c0), notes: 'Constant learning rate and growth; single-product plant; costs in constant money.' },
  };
}

// ---------- mode: material passport and end of life ----------
function fnv(s) { let h = 0x811c9dc5; for (let i = 0; i < s.length; i++) { h ^= s.charCodeAt(i); h = Math.imul(h, 0x01000193) >>> 0; } return ('0000000' + h.toString(16)).slice(-8); }
function runPassport(inp, ctx) {
  const w = [], bus = (ctx && ctx.bus) || {}, rt = ROUTES[inp.route] || ROUTES.zno, txt = v => typeof v === 'string' && v.trim() ? v.trim() : null;
  let formula = txt(inp.pFormula) || (typeof bus['material.formula'] === 'string' ? bus['material.formula'] : null) || rt.pr[0].f, comp = {}, M = 0;
  try { comp = chem.massFractions(formula); M = chem.molarMass(formula); } catch (e) { w.push(`Formula "${formula}" could not be parsed (${e.message}); the product of the selected route is used.`); formula = rt.pr[0].f; comp = chem.massFractions(formula); M = chem.molarMass(formula); }
  const name = txt(inp.pName) || (typeof bus['material.name'] === 'string' ? bus['material.name'] : null) || rt.name, pick = pre => Object.fromEntries(Object.entries(bus).filter(([k, v]) => k.startsWith(pre) && (typeof v === 'number' ? Number.isFinite(v) : typeof v === 'string')).map(([k, v]) => [k.slice(pre.length), v]));
  const props = { ...pick('prop.'), ...Object.fromEntries(Object.entries(pick('eff.')).map(([k, v]) => ['effective_' + k, v])) }, geom = pick('geom.'), proc = pick('proc.'), cred = pick('cred.'), risk = pick('risk.');
  const els = Object.keys(comp), extra = []; String(inp.elements || '').split(/[;,\s]+/).filter(Boolean).forEach(f => { try { Object.keys(chem.parseFormula(f).atoms).forEach(e => extra.push(e)); } catch { w.push(`Additional formula "${f}" could not be read.`); } });
  const crm = [...new Set([...els, ...extra])].filter(e => CRM[e]), co2 = Number.isFinite(bus['econ.co2_per_kg']) ? bus['econ.co2_per_kg'] : inp.co2Primary;
  // end-of-life scenarios per kg of end-of-life product
  const g = inp.grid / 1e3, wC = comp.C || 0, metal = ['Au', 'Ag', 'Cu', 'Zn'].find(e => comp[e]), mv = metal ? metalPrice(inp, metal) * comp[metal] : 0, val = inp.valRec > 0 ? inp.valRec : mv, rec = inp.recOut / 100 * inp.recEff / 100, recProc = inp.recEff / 100;
  const eol = {
    landfill: { label: 'Landfill', cost: inp.feeLand / 1e3, co2: 0.02 + 0.05 * wC * 44 / 12, rec: 0 },
    incineration: { label: 'Incineration with energy recovery', cost: inp.feeInc / 1e3, co2: wC * 44 / 12 - inp.lhv / 3.6 * 0.25 * g, rec: 0 },
    recycling: { label: 'Material recycling', cost: inp.costRec - recProc * val, co2: inp.eRec * g - recProc * co2, rec: recProc },
  }, order = Object.entries(eol), bestC = order.slice().sort((a, b) => a[1].cost - b[1].cost)[0], bestE = order.slice().sort((a, b) => a[1].co2 - b[1].co2)[0], decl = eol[inp.eolRoute] || eol.recycling;
  if (inp.valRec <= 0 && metal) w.push(`Recovered-material value taken from the ${metal} market price and the ${pct(comp[metal])} metal content (${F(mv, 3)} USD/kg); enter a value to override.`);
  if (decl !== bestE[1]) w.push(`The declared end-of-life route (${decl.label.toLowerCase()}) is not the route with the lowest greenhouse-gas emissions in this comparison (${bestE[1].label.toLowerCase()}).`);
  w.push('The passport assembles values supplied by the other suites and by the user; it records them with their provenance but does not verify them. Fields reported as missing must be completed from qualified sources, not estimated.');
  // passport document
  const doc = {
    schema: 'designer-materials/material-passport/1.0',
    identity: { name, formula, formulaDisplay: chem.prettyFormula(formula), molarMass_g_mol: +M.toFixed(4), materialClass: typeof bus['material.class'] === 'string' ? bus['material.class'] : null, productionRoute: rt.label },
    composition: { massFractions: Object.fromEntries(Object.entries(comp).map(([e, v]) => [e, +v.toFixed(5)])), additionalElements: [...new Set(extra)], criticalRawMaterialScreening: crm.map(e => ({ element: e, listedAs: CRM[e] })) },
    processRecipe: { route: rt.label, reaction: rt.physical ? null : chem.ionicEquations(rt.re, rt.pr).molecular, parameters: proc },
    properties: { intrinsic: props, geometry: geom },
    hazards: { classification: txt(inp.hazText), exposure: risk, note: 'Hazard classification must be taken from the safety data sheet; modelling does not establish safety.' },
    environment: { carbonFootprint_kgCO2e_per_kg: Number.isFinite(co2) ? co2 : null, footprintSource: Number.isFinite(bus['econ.co2_per_kg']) ? 'life-cycle screening (Suite 28)' : 'user entry', gridIntensity_gCO2_per_kWh: inp.grid },
    endOfLife: { declaredRoute: inp.eolRoute, recoveryRate: +(inp.eolRoute === 'recycling' ? rec : 0).toFixed(4), scenarios: Object.fromEntries(order.map(([k, v]) => [k, { cost_USD_per_kg: +v.cost.toFixed(4), co2_kg_per_kg: +v.co2.toFixed(4) }])) },
    credibility: { class: typeof cred.class === 'string' ? cred.class : null, evidence: cred, note: 'Model-derived property values carry the credibility class of the model that produced them.' },
    provenance: { batch: txt(inp.batch), producer: txt(inp.producer), issued: txt(inp.issued), generator: 'Designer Materials platform, Suite 28' },
  };
  const fields = [['Name', doc.identity.name], ['Formula', doc.identity.formula], ['Material class', doc.identity.materialClass], ['Composition', els.length ? 1 : null], ['Process parameters', Object.keys(proc).length ? 1 : null], ['Intrinsic or effective properties (≥ 3)', Object.keys(props).length >= 3 ? 1 : null], ['Geometry or morphology', Object.keys(geom).length ? 1 : null], ['Hazard classification', doc.hazards.classification], ['Exposure or risk assessment', Object.keys(risk).length ? 1 : null], ['Carbon footprint', Number.isFinite(co2) && co2 > 0 ? 1 : null], ['End-of-life route', inp.eolRoute], ['Credibility class of model data', doc.credibility.class], ['Batch identifier', doc.provenance.batch], ['Producer', doc.provenance.producer], ['Date of issue', doc.provenance.issued]];
  const filled = fields.filter(f => f[1] !== null && f[1] !== undefined).length, compl = filled / fields.length, missing = fields.filter(f => f[1] === null || f[1] === undefined).map(f => f[0]), nBus = Object.keys(props).length + Object.keys(geom).length + Object.keys(proc).length + Object.keys(cred).length + Object.keys(risk).length;
  const body = JSON.stringify(doc), sum = fnv(body); doc.integrity = { checksum_fnv1a32: sum, completeness: +compl.toFixed(3), missingFields: missing };
  const rr = N.linspace(0, 1, 21), eqs = routeEquations(inventory({ ...inp, yieldFrac: 0.9, excess: 0, cPrec: 50, washVol: 0, solvType: 'water', solvRec: 0, T: 20, heatLoss: 0, heatRec: 0, eExtra: 0 }));
  eqs.push({ kind: 'formula', label: 'Composition', text: `${chem.prettyFormula(formula)}:  ` + Object.entries(comp).map(([e, v]) => `w(${e}) = ${F(100 * v, 4)} %`).join(', '), note: `Molar mass ${F(M, 6)} g/mol.` },
    { kind: 'formula', label: 'End-of-life cost and emissions per kg', text: 'recycling: c = c_process − η·v_recovered,  e = E_rec·EF_grid − η·GWP_primary;   incineration: e = w_C·(44/12) − LHV·η_el·EF_grid;   landfill: c = gate fee' });
  if (wC > 0) eqs.push({ kind: 'chemical', label: 'Carbon oxidation on incineration', text: chem.equationString(['C', 'O2'], ['CO2']), note: `Carbon mass fraction ${pct(wC)} gives ${F(wC * 44 / 12, 3)} kg CO₂ per kg.` });
  return {
    summary: `The passport for ${name} (${chem.prettyFormula(formula)}, ${F(M, 5)} g/mol) contains ${filled} of ${fields.length} required fields (${pct(compl)} complete) and ${nBus} values taken from other suites; checksum ${sum}. ${missing.length ? `Missing: ${missing.slice(0, 5).join(', ')}${missing.length > 5 ? '…' : ''}.` : 'No required field is missing.'} Per kilogram at end of life, recycling costs ${F(eol.recycling.cost, 3)} USD and emits ${F(eol.recycling.co2, 3)} kg CO₂e (negative values are net credits), incineration ${F(eol.incineration.cost, 3)} USD and ${F(eol.incineration.co2, 3)} kg, landfill ${F(eol.landfill.cost, 3)} USD and ${F(eol.landfill.co2, 3)} kg; the lowest-emission route is ${bestE[1].label.toLowerCase()}.`,
    kpis: [
      kp('compl', 'Passport completeness', 100 * compl, '%', `Share of the ${fields.length} required passport fields that hold a value from the user or from other suites; missing fields are listed in the passport and must not be filled with estimates.`, compl >= 0.8),
      kp('nbus', 'Values taken from other suites', nBus, '', 'Number of property, geometry, process, risk and credibility values read from the shared data bus; run the relevant suites first to populate them.', nBus > 0),
      kp('checksum', 'Integrity checksum', sum, '', 'FNV-1a 32-bit checksum of the passport content: any later change of a recorded value changes the checksum. It detects accidental alteration and is not a cryptographic signature.', null),
      kp('crm', 'Critical raw materials flagged', crm.length, '', crm.length ? `Elements on the screening list: ${crm.join(', ')}; recorded in the composition section for recovery planning.` : 'No element of the declared composition is on the built-in critical-raw-material screening list.', crm.length === 0),
      kp('c_rec', 'Recycling: net cost', eol.recycling.cost, 'USD/kg', `Process cost of ${inp.costRec} USD/kg less the value of recovered material (${F(val, 3)} USD/kg at ${inp.recEff} % process efficiency); negative values mean recycling pays for itself.`, eol.recycling.cost <= eol.landfill.cost),
      kp('e_rec', 'Recycling: net emissions', eol.recycling.co2, 'kg CO₂e/kg', `Recycling energy (${inp.eRec} kWh/kg at ${inp.grid} g/kWh) less the avoided primary production (${F(co2, 3)} kg CO₂e/kg × recovery); negative values are a net credit.`, eol.recycling.co2 < 0),
      kp('c_inc', 'Incineration: cost', eol.incineration.cost, 'USD/kg', 'Gate fee for incineration with energy recovery, per kilogram.', null),
      kp('e_inc', 'Incineration: net emissions', eol.incineration.co2, 'kg CO₂e/kg', `CO₂ from the carbon content (${pct(wC)}) less the credit for electricity generated at 25 % efficiency from the heating value; inorganic residues go to ash.`, eol.incineration.co2 <= 0.05),
      kp('c_land', 'Landfill: cost', eol.landfill.cost, 'USD/kg', 'Landfill gate fee per kilogram; it carries no recovery of material or energy.', null),
      kp('e_land', 'Landfill: emissions', eol.landfill.co2, 'kg CO₂e/kg', 'Transport and site operation (0.02 kg) plus slow degradation of a small part of any organic carbon; the material value and embodied emissions are lost.', null),
      kp('best', 'Lowest-emission end-of-life route', bestE[1].label, '', `Route with the lowest net greenhouse-gas emissions; the cheapest route is ${bestC[1].label.toLowerCase()}. The declared route is ${decl.label.toLowerCase()}.`, decl === bestE[1]),
      kp('rec_rate', 'Overall material recovery rate', rec, '–', `Collection (${inp.recOut} %) × recycling efficiency (${inp.recEff} %): the share of the material placed on the market that returns as secondary raw material.`, rec > 0.5),
    ],
    series: [
      { key: 'eol', title: 'End-of-life routes compared', xlabel: 'route (1 landfill, 2 incineration, 3 recycling)', ylabel: 'per kg of product', lines: [{ name: 'cost (USD/kg)', x: [1, 2, 3], y: order.map(o => o[1].cost) }, { name: 'net emissions (kg CO₂e/kg)', x: [1, 2, 3], y: order.map(o => o[1].co2) }], explain: 'Cost and net greenhouse-gas emissions of the three end-of-life routes per kilogram of product. Negative values are credits from recovered material or energy.' },
      { key: 'rr', title: 'Recycling balance versus process efficiency', xlabel: 'recycling process efficiency (–)', ylabel: 'per kg of product', lines: [{ name: 'net emissions (kg CO₂e/kg)', x: rr, y: rr.map(e => inp.eRec * g - e * co2) }, { name: 'net cost (USD/kg)', x: rr, y: rr.map(e => inp.costRec - e * val) }, { name: 'zero', x: [0, 1], y: [0, 0] }], explain: 'Net emissions and net cost of recycling as the recovered fraction rises. The crossing of zero is the minimum efficiency at which recycling gives an environmental or economic benefit.' },
    ],
    tables: [
      { key: 'fields', title: 'Passport fields', columns: ['Field', 'Status'], rows: fields.map(f => [f[0], f[1] !== null && f[1] !== undefined ? 'present' : 'missing']), explain: 'Required fields of the digital material passport and whether each holds a value. Fields fed by other suites become present once those suites have been run.' },
      { key: 'comp', title: 'Elemental composition', columns: ['Element', 'Mass fraction (%)', 'Critical-material screening'], rows: Object.entries(comp).map(([e, v]) => [e, F(100 * v, 4), CRM[e] || 'not listed']), explain: 'Mass fractions from the chemical formula with the critical-raw-material screening result for each element.' },
      { key: 'eol', title: 'End-of-life scenarios', columns: ['Route', 'Cost (USD/kg)', 'Net emissions (kg CO₂e/kg)', 'Material recovered (kg/kg)'], rows: order.map(([, v]) => [v.label, F(v.cost, 4), F(v.co2, 4), F(v.rec, 3)]), explain: 'Screening comparison of landfill, incineration and recycling. Gate fees, energy and credits are user-adjustable; collection losses apply before any of the routes.' },
    ],
    equations: eqs,
    checks: [
      { name: 'Composition sums to 100 %', pass: Math.abs(N.sum(Object.values(comp)) - 1) < 1e-9, detail: `Σ w = ${F(N.sum(Object.values(comp)), 8)}` },
      { name: 'Passport at least 80 % complete', pass: compl >= 0.8, detail: `${filled} of ${fields.length} fields` },
      { name: 'Provenance recorded (batch, producer, date)', pass: !!(doc.provenance.batch && doc.provenance.producer && doc.provenance.issued), detail: [doc.provenance.batch, doc.provenance.producer, doc.provenance.issued].map(v => v || 'missing').join(' / ') },
      { name: 'Hazard classification recorded', pass: !!doc.hazards.classification, detail: doc.hazards.classification || 'missing: take from the safety data sheet' },
      { name: 'Declared end-of-life route has the lowest emissions', pass: decl === bestE[1], detail: `declared ${decl.label.toLowerCase()}, lowest ${bestE[1].label.toLowerCase()}` },
      { name: 'Recycling gives a net emission credit', pass: eol.recycling.co2 < 0, detail: `${F(eol.recycling.co2, 4)} kg CO₂e/kg` },
    ],
    warnings: w,
    recommendations: [
      { priority: missing.length ? 'high' : 'low', text: missing.length ? `Complete the missing passport fields (${missing.join(', ')}): run the property, safety, validation and life-cycle suites for model-derived values and take hazard classification and provenance from controlled documents.` : 'The passport is complete; issue it with the batch record and re-issue it whenever a recorded value or the recipe changes.', basis: 'Passport completeness check.' },
      { priority: 'medium', text: eol.recycling.co2 < Math.min(eol.landfill.co2, eol.incineration.co2) ? `Declare recycling as the end-of-life route and design for it: at ${inp.recEff} % process efficiency each kilogram recycled avoids ${F(-eol.recycling.co2, 3)} kg CO₂e, but only ${pct(rec)} of the material returns unless collection rises above ${inp.recOut} %.` : `Recycling does not pay back environmentally with these figures (process energy ${inp.eRec} kWh/kg against a primary footprint of ${F(co2, 3)} kg CO₂e/kg); reduce the recycling energy or reserve recycling for the valuable fraction.`, basis: 'End-of-life scenario comparison.', sustainability: true },
      { priority: crm.length ? 'high' : 'low', text: crm.length ? `Record the mass and location of ${crm.join(', ')} in the product so that recyclers can target them, and include a take-back clause for products that contain them.` : 'Keep the composition record current: adding dopants or electrodes can introduce critical raw materials that change the recovery strategy.', basis: 'Critical-raw-material screening of the composition.', sustainability: true },
      { priority: 'medium', text: 'Attach the credibility class from Suite 24 to every model-derived property in the passport and distinguish measured from predicted values; recipients must be able to see which figures rest on validated evidence.', basis: 'Provenance and credibility fields of the passport.' },
    ],
    publish: { 'passport.completeness': compl, 'passport.checksum': sum, 'econ.eol_cost_per_kg': decl.cost, 'econ.eol_co2_per_kg': decl.co2, 'material.formula': formula, 'material.name': name },
    exports: [{ name: 'material_passport.json', mime: 'application/json', content: JSON.stringify(doc, null, 1), description: 'Machine-readable digital material passport: identity, composition, process recipe, properties, hazards, carbon footprint, end-of-life route, credibility and provenance, with integrity checksum.' }],
    numerics: { method: 'Assembly of bus and user values into a structured record; composition from the chem kernel; algebraic end-of-life cost and emission balances; FNV-1a checksum', dof: fields.length, steps: 1, iterations: 1, residual: Math.abs(N.sum(Object.values(comp)) - 1), notes: 'End-of-life figures are per kilogram of product reaching the route; avoided-burden credits follow the substitution approach.' },
  };
}

export default {
  id: 28,
  tagline: 'Estimates production cost, investment return, life-cycle impacts, circularity and learning effects of a nanomaterial or metamaterial process, and assembles its digital material passport.',
  inputs,
  presets: [
    { name: 'ZnO nanorods, hydrothermal production 20 t/yr', note: 'Dilute aqueous route: precursor cost indexed to the zinc price, solvent and heating energy dominate the inventory.', values: { mode: 'tea', route: 'zno', cap: 20000, capRef: 20000, cPrec: 200, T: 90, premium: 4, pReag: 0.5, price: 60, solvType: 'water', solvRec: 50 } },
    { name: 'Ag nanoparticle colloid, 200 kg/yr', note: 'Noble-metal route where the silver price dominates the unit cost and the footprint; small plant, high-value product.', values: { mode: 'tea', route: 'ag', cap: 200, capRef: 200, eqRef: 150000, cPrec: 5, T: 25, premium: 1.3, pReag: 30, price: 4000, excess: 50, nOp: 1, opHours: 2000, loading: 0.05, eExtra: 20 } },
    { name: 'Au nanoparticle colloid: life-cycle screening', note: 'Gold precursor dominates the footprint; functional unit of one sensor device with microgram loading.', values: { mode: 'lca', route: 'au', cPrec: 1, T: 100, efPrec: 12000, efReag: 2, cedPrec: 180000, cedReag: 30, fu: 'device', loading: 0.02, devArea: 0.25, recIn: 30, recOut: 60, recEff: 95, yieldFrac: 0.95, eSave: 0 } },
    { name: 'TiO₂ metasurface layer by sol–gel: cost per m²', note: 'Alkoxide route in alcohol with solvent recovery; cost expressed per square metre of patterned metasurface and per device.', values: { mode: 'tea', route: 'tio2', cap: 2000, capRef: 2000, eqRef: 250000, useIndex: false, pPrec: 8, pReag: 0, solvType: 'isopropanol', pSolv: 1200, solvRec: 85, cPrec: 300, T: 60, price: 400, loading: 1.5, devArea: 4, integ: 2.5, eExtra: 15 } },
    { name: '3-D printed auxetic lattice: investment case', note: 'Photopolymer additive-manufacturing line for architected lattice parts: discounted cash flow with ramp-up and uncertainty.', values: { mode: 'dcf', capex: 9e5, opex: 3.2e5, cap: 6000, price: 110, varFrac: 55, ramp: 3, life: 8, rate: 12, cvVol: 25, cvPrice: 15 } },
    { name: '3-D printed lattice: passport and end of life', note: 'Acrylate photopolymer lattice: carbon-containing product where incineration emits CO₂ and recycling is limited.', values: { mode: 'passport', route: 'resin', lhv: 27, costRec: 1.2, eRec: 1.5, co2Primary: 6, recOut: 20, recEff: 60, eolRoute: 'incineration' } },
    { name: 'ZnO nanorods: learning curve to commodity cost', note: 'Experience-curve projection from pilot volumes, with scale-up and make-or-buy for a coating manufacturer.', values: { mode: 'learning', c0: 45, q0: 5000, lr: 18, cTarget: 15, qT: 5e5, qYear: 3000, growth: 50 } },
  ],
  run(inp, ctx) {
    const r = ({ tea: runTea, dcf: runDcf, lca: runLca, learning: runLearning, passport: runPassport }[inp.mode] || runTea)(inp, ctx);
    for (const k of r.kpis) if (typeof k.value === 'number' && !Number.isFinite(k.value)) { r.warnings.push(`${k.label} could not be evaluated for these inputs and is reported as 0.`); k.value = 0; k.good = null; }
    r.series = (r.series || []).map(s => ({ ...s, lines: s.lines.map(l => { const x = [], y = []; for (let i = 0; i < l.x.length; i++) if (Number.isFinite(l.x[i]) && Number.isFinite(l.y[i]) && (!s.logy || l.y[i] > 0) && (!s.logx || l.x[i] > 0)) { x.push(l.x[i]); y.push(l.y[i]); } return { ...l, x, y }; }).filter(l => l.x.length) })).filter(s => s.lines.length);
    for (const k of Object.keys(r.publish)) if (typeof r.publish[k] === 'number' && !Number.isFinite(r.publish[k])) delete r.publish[k];
    r.warnings.push('Screening-level economic and environmental estimates from user-adjustable prices and factors: use them to compare options and find hot spots, and replace the factors by quotations and qualified data before any commitment or external claim.');
    return r;
  },
};
