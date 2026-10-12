// Reference benchmarks for Suite 28 — economics, sustainability, data and lifecycle.
// Expected values are closed-form identities of engineering economics, stoichiometry from standard atomic
// weights, and hand-worked cases; none of them is taken from the suite's own output.

// --- standard atomic weights (IUPAC abridged)
const AW = { H: 1.008, C: 12.011, N: 14.007, O: 15.999, Na: 22.990, Ti: 47.867, Zn: 65.38, Ag: 107.868 };
const M_ZNO = AW.Zn + AW.O, M_NO3 = AW.N + 3 * AW.O, M_ZN_NIT = AW.Zn + 2 * M_NO3, M_H2O = 2 * AW.H + AW.O;
const M_ZN_NIT_6W = M_ZN_NIT + 6 * M_H2O, M_NAOH = AW.Na + AW.O + AW.H, M_NANO3 = AW.Na + M_NO3;
const M_AGNO3 = AW.Ag + M_NO3, M_TIO2 = AW.Ti + 2 * AW.O, M_TTIP = AW.Ti + 4 * (AW.O + 3 * AW.C + 7 * AW.H);
const M_MMA = 5 * AW.C + 8 * AW.H + 2 * AW.O, M_CO2 = AW.C + 2 * AW.O;

// --- engineering economics
const annuity = (i, n) => (1 - Math.pow(1 + i, -n)) / i;            // present worth of a unit annuity
const crf = (i, n) => 1 / annuity(i, n);                             // capital recovery factor
// investment case without tax, working capital or ramp-up: 1.0 M USD, 10 000 kg/yr at 50 USD/kg, operating cost 0.2 M USD/yr
const DCF = { mode: 'dcf', capex: 1e6, cap: 10000, price: 50, opex: 2e5, varFrac: 45, ramp: 0, tax: 0, depYears: 5, wc: 0, life: 10, rate: 10, cvPrice: 0, cvVol: 0, cvOpex: 0, cvCapex: 0, nMc: 200, seed: 1 };
const REV = 5e5, OP = 2e5, A10 = annuity(0.1, 10), NPV0 = -1e6 + (REV - OP) * A10;
// lognormal price multiplier with mean one: x = exp(s·z − s²/2), s² = ln(1 + CV²)
const erf = x => { const t = 1 / (1 + 0.3275911 * Math.abs(x)), y = 1 - (((((1.061405429 * t - 1.453152027) * t) + 1.421413741) * t - 0.284496736) * t + 0.254829592) * t * Math.exp(-x * x); return x >= 0 ? y : -y; };
const Phi = z => 0.5 * (1 + erf(z / Math.SQRT2));
const S30 = Math.sqrt(Math.log(1 + 0.3 * 0.3)), M_STAR = (1e6 / A10 + OP) / REV;   // price multiplier at which NPV = 0

// --- hand-worked techno-economic case: photopolymer route, feedstock only
const TEA_HAND = { mode: 'tea', route: 'resin', yieldFrac: 1, washVol: 0, eExtra: 0, useIndex: false, pPrec: 10, nOp: 0, maint: 0, cap: 20000, capRef: 20000, eqRef: 4e5, lang: 'mixed', solvRec: 50, rate: 10, life: 10, price: 60 };
const CAPEX_REF = 3.63 * 4e5, UNIT_HAND = 10 + crf(0.1, 10) * CAPEX_REF / 20000;

// --- hand-worked ZnO inventory: stoichiometric, concentrated, no washing, no heating
const INV = { route: 'zno', yieldFrac: 1, excess: 0, cPrec: 1000, washVol: 0, solvType: 'water', solvRec: 0, T: 20, heatLoss: 0, heatRec: 0, eExtra: 3 };
const N_ZNO = 1000 / M_ZNO;                                          // mol per kg of product
const m_prec = M_ZN_NIT_6W / M_ZNO, m_naoh = 2 * M_NAOH / M_ZNO, m_salt = 2 * M_NANO3 / M_ZNO, m_water = N_ZNO / 1000 * 1000;   // kg per kg ZnO
const GWP_HAND = m_prec * 3 + m_naoh * 1.2 + m_water * 0.0005 + 3 * 0.35 + m_salt * 0.1 + (m_prec + m_naoh) / 1000 * 500 * 0.1;
const LCA_HAND = { mode: 'lca', ...INV, fu: 'kg', efPrec: 3, efReag: 1.2, efSolv: 0.0005, grid: 350, efWaste: 0.1, transKm: 500, efTrans: 0.1 };

// --- material circularity indicator, hand-worked: 50 % recycled feedstock, 50 % collected, 80 % recycling efficiency
const MCI = (() => { const FR = 0.5, CR = 0.5, E = 0.8, V = 1 - FR, WC = (1 - E) * CR, WF = (1 - E) * FR / E, W = 1 - CR + (WF + WC) / 2, LFI = (V + W) / (2 + (WF - WC) / 2); return { LFI, MCI: 1 - 0.9 * LFI }; })();   // LFI = 1.1125/2.0125

export default [
  // ------------------------------------------------------------------ techno-economic analysis
  { id: 's28-zno-precursor-demand', suite: 28, title: 'Zinc nitrate hexahydrate demand per kg ZnO at 90 % yield', mode: 'tea',
    inputs: { mode: 'tea', route: 'zno', yieldFrac: 0.9 }, kpi: 'm_prec', expected: M_ZN_NIT_6W / M_ZNO / 0.9, tol: 5e-4, kind: 'analytical',
    reference: 'Stoichiometry of Zn(NO₃)₂·6H₂O + 2 NaOH → ZnO + 2 NaNO₃ + 7 H₂O with IUPAC standard atomic weights', note: '297.48/81.38/0.9 = 4.06 kg/kg.' },
  { id: 's28-zno-atom-economy', suite: 28, title: 'Atom economy of ZnO precipitation from zinc nitrate and sodium hydroxide', mode: 'tea',
    inputs: { mode: 'tea', route: 'zno' }, kpi: 'atom', expected: M_ZNO / (M_ZN_NIT + 2 * M_NAOH), tol: 5e-4, kind: 'analytical',
    reference: 'Trost atom economy, M_product/ΣM_reactants for Zn(NO₃)₂ + 2 NaOH → ZnO + 2 NaNO₃ + H₂O', note: '81.38/269.38 = 0.302.' },
  { id: 's28-zno-e-factor', suite: 28, title: 'E-factor of stoichiometric ZnO precipitation (sodium nitrate by-product)', mode: 'tea',
    inputs: { mode: 'tea', ...INV }, kpi: 'efactor', expected: m_salt, tol: 5e-4, kind: 'analytical',
    reference: 'Sheldon E-factor (water excluded): 2 mol NaNO₃ per mol ZnO at full yield and no reagent excess', note: '2 × 84.99/81.38 = 2.09 kg/kg.' },
  { id: 's28-ag-precursor-demand', suite: 28, title: 'Silver nitrate demand per kg of silver at full yield', mode: 'tea',
    inputs: { mode: 'tea', route: 'ag', yieldFrac: 1 }, kpi: 'm_prec', expected: M_AGNO3 / AW.Ag, tol: 5e-4, kind: 'analytical',
    reference: 'Stoichiometry: one mole of AgNO₃ (169.87 g/mol) per mole of Ag (107.87 g/mol)' },
  { id: 's28-tio2-atom-economy', suite: 28, title: 'Atom economy of TiO₂ from titanium isopropoxide hydrolysis', mode: 'tea',
    inputs: { mode: 'tea', route: 'tio2', useIndex: false }, kpi: 'atom', expected: M_TIO2 / (M_TTIP + 2 * M_H2O), tol: 5e-4, kind: 'analytical',
    reference: 'Stoichiometry of Ti(OC₃H₇)₄ + 2 H₂O → TiO₂ + 4 C₃H₇OH', note: '79.87/320.25 = 0.249.' },
  { id: 's28-heating-energy', suite: 28, title: 'Sensible heat to bring the reaction liquid from 20 °C to 90 °C', mode: 'tea',
    inputs: { mode: 'tea', ...INV, T: 90, eExtra: 0 }, kpi: 'energy', expected: m_water * 4180 * 70 / 3.6e6, tol: 1e-3, kind: 'analytical',
    reference: 'Q = m c_p ΔT for water (c_p = 4.18 kJ/kg/K, 1000 kg/m³) at 1 mol/L of product', note: '12.29 kg of water per kg ZnO: 1.0 kWh/kg.' },
  { id: 's28-six-tenths-rule', suite: 28, title: 'Six-tenths rule: capital cost of a plant ten times the reference capacity', mode: 'tea',
    inputs: { mode: 'tea', cap: 200000, capRef: 20000, eqRef: 4e5, scaleExp: 0.6, lang: 'mixed', solvRec: 50 }, kpi: 'capex', expected: CAPEX_REF * Math.pow(10, 0.6), tol: 1e-6, kind: 'analytical',
    reference: 'Cost–capacity relation C = C_ref(S/S_ref)^0.6 and Lang factor 3.63 for solid–fluid plants (Peters and Timmerhaus, Plant Design and Economics for Chemical Engineers)', note: '10^0.6 = 3.981.' },
  { id: 's28-lang-fluid', suite: 28, title: 'Lang factor for a fluid-processing plant', mode: 'tea',
    inputs: { mode: 'tea', cap: 20000, capRef: 20000, eqRef: 4e5, lang: 'fluid', solvRec: 50 }, kpi: 'capex', expected: 4.74 * 4e5, tol: 1e-6, kind: 'reference-data',
    reference: 'Lang factors for fixed capital investment: 3.10 solid, 3.63 solid–fluid, 4.74 fluid processing (Peters and Timmerhaus)' },
  { id: 's28-levelised-unit-cost', suite: 28, title: 'Levelised unit cost: feedstock plus annuity capital charge (hand-worked)', mode: 'tea',
    inputs: TEA_HAND, kpi: 'unit', expected: UNIT_HAND, tol: 1e-6, kind: 'analytical',
    reference: 'Capital recovery factor CRF = i(1 + i)ⁿ/((1 + i)ⁿ − 1) (Peters and Timmerhaus)', note: 'CRF(10 %, 10 yr) = 0.16275 on 1.452 M USD over 20 t/yr adds 11.82 USD/kg to 10 USD/kg of feedstock.' },
  { id: 's28-break-even-volume', suite: 28, title: 'Break-even volume of the hand-worked case', mode: 'tea',
    inputs: TEA_HAND, kpi: 'be', expected: crf(0.1, 10) * CAPEX_REF / (60 - 10), tol: 1e-6, kind: 'analytical',
    reference: 'Break-even analysis Q_BE = fixed cost/(price − variable cost)' },

  // ------------------------------------------------------------------ discounted cash flow
  { id: 's28-npv-annuity', suite: 28, title: 'NPV of a level cash flow equals the annuity closed form', mode: 'dcf',
    inputs: DCF, kpi: 'npv', expected: NPV0, tol: 1e-9, kind: 'analytical',
    reference: 'Present worth of a uniform series, P = A[1 − (1 + i)⁻ⁿ]/i (Brealey and Myers, Principles of Corporate Finance)', note: '0.3 M USD/yr for 10 years at 10 %: NPV = 843 370 USD.' },
  { id: 's28-npv-zero-rate', suite: 28, title: 'NPV at zero discount rate equals the undiscounted sum', mode: 'dcf',
    inputs: { ...DCF, rate: 0 }, kpi: 'npv', expected: -1e6 + 10 * (REV - OP), tol: 1e-9, kind: 'conservation',
    reference: 'Definition of net present value at i = 0' },
  { id: 's28-irr-identity', suite: 28, title: 'IRR of an investment priced at the 20 % annuity factor', mode: 'dcf',
    inputs: { ...DCF, life: 5, capex: (REV - OP) * annuity(0.2, 5) }, kpi: 'irr', expected: 20, tol: 1e-6, kind: 'analytical',
    reference: 'Internal rate of return: NPV(IRR) = 0 (Brealey and Myers, Principles of Corporate Finance)', note: 'The investment is set to 0.3 M USD × a(20 %, 5 yr) = 897 184 USD, so the IRR is 20 % by construction.' },
  { id: 's28-zero-npv-price', suite: 28, title: 'Selling price for zero NPV equals operating cost plus capital recovery', mode: 'dcf',
    inputs: DCF, kpi: 'p_be', expected: (OP + 1e6 * crf(0.1, 10)) / 10000, tol: 1e-6, kind: 'analytical',
    reference: 'Capital recovery factor: annual revenue requirement = OPEX + CRF·CAPEX' },
  { id: 's28-profitability-index', suite: 28, title: 'Profitability index of the level-cash-flow case', mode: 'dcf',
    inputs: DCF, kpi: 'pi', expected: (REV - OP) * A10 / 1e6, tol: 1e-9, kind: 'analytical',
    reference: 'Profitability index = present value of inflows/investment' },
  { id: 's28-equivalent-annual-value', suite: 28, title: 'Equivalent annual value: cash flow less the capital recovery charge', mode: 'dcf',
    inputs: DCF, kpi: 'annuity', expected: (REV - OP) - 1e6 * crf(0.1, 10), tol: 1e-9, kind: 'analytical',
    reference: 'Equivalent annual worth, NPV × CRF' },
  { id: 's28-npv-tax-depreciation', suite: 28, title: 'After-tax NPV with straight-line depreciation (hand-worked)', mode: 'dcf',
    inputs: { ...DCF, life: 5, depYears: 5, tax: 25 }, kpi: 'npv', expected: -1e6 + ((REV - OP - 2e5) * 0.75 + 2e5) * annuity(0.1, 5), tol: 1e-9, kind: 'analytical',
    reference: 'After-tax cash flow CF = (R − O − D)(1 − τ) + D (Brealey and Myers, Principles of Corporate Finance)', note: 'CF = 275 000 USD/yr for 5 years: NPV = 42 466 USD.' },
  { id: 's28-npv-working-capital', suite: 28, title: 'NPV with working capital invested in year 1 and recovered in the last year', mode: 'dcf',
    inputs: { ...DCF, wc: 10 }, kpi: 'npv', expected: NPV0 - 0.1 * REV / 1.1 + 0.1 * REV / Math.pow(1.1, 10), tol: 1e-9, kind: 'analytical',
    reference: 'Discounting of working-capital outlay and recovery (Brealey and Myers, Principles of Corporate Finance)' },
  { id: 's28-monte-carlo-degenerate', suite: 28, title: 'Monte Carlo with zero uncertainty collapses to the deterministic NPV', mode: 'dcf',
    inputs: DCF, kpi: 'mc_p50', expected: NPV0, tol: 1e-9, kind: 'conservation',
    reference: 'Consistency of the sampling with the base case when all coefficients of variation are zero' },
  { id: 's28-monte-carlo-median', suite: 28, title: 'Monte Carlo median NPV under lognormal price uncertainty (CV 30 %)', mode: 'dcf',
    inputs: { ...DCF, cvPrice: 30, nMc: 30000 }, kpi: 'mc_p50', expected: -1e6 + (REV * Math.exp(-S30 * S30 / 2) - OP) * A10, tol: 0.03, kind: 'analytical',
    reference: 'Median of a mean-one lognormal variable, exp(−σ²/2) with σ² = ln(1 + CV²); NPV is linear in price without tax',
    note: 'Statistical benchmark with 30 000 seeded samples: the standard error of the sample median is about 0.9 % of the expected value, so 3 % is roughly three standard errors.' },
  { id: 's28-monte-carlo-loss-probability', suite: 28, title: 'Monte Carlo probability of negative NPV under lognormal price uncertainty (CV 30 %)', mode: 'dcf',
    inputs: { ...DCF, cvPrice: 30, nMc: 30000 }, kpi: 'p_neg', expected: Phi((Math.log(M_STAR) + S30 * S30 / 2) / S30), abs: 0.008, kind: 'analytical',
    reference: 'Lognormal distribution function: P(x < x*) = Φ[(ln x* + σ²/2)/σ] for a mean-one multiplier',
    note: 'Exact value about 0.172; the binomial standard error for 30 000 samples is 0.0022, so the absolute tolerance is about 3.6 standard errors.' },

  // ------------------------------------------------------------------ life-cycle screening and circularity
  { id: 's28-gwp-hand-inventory', suite: 28, title: 'Cradle-to-gate footprint of a hand-worked ZnO inventory', mode: 'lca',
    inputs: LCA_HAND, kpi: 'gwp', expected: GWP_HAND, tol: 1e-3, kind: 'analytical',
    reference: 'Process-based inventory, GWP = Σ mᵢ·EFᵢ + E·EF_grid + m_waste·EF_waste + (t·km)·EF_transport, with stoichiometric masses from standard atomic weights',
    note: 'Per kg ZnO: 3.655 kg nitrate hexahydrate, 0.983 kg NaOH, 12.29 kg water, 3 kWh, 2.089 kg salt waste, 2.32 t·km: 13.64 kg CO₂e with the stated screening factors.' },
  { id: 's28-gwp-per-square-metre', suite: 28, title: 'Functional-unit conversion: footprint per m² at 2 g/m² loading', mode: 'lca',
    inputs: { ...LCA_HAND, fu: 'm2', loading: 2 }, kpi: 'gwp', expected: GWP_HAND * 2e-3, tol: 1e-3, kind: 'analytical',
    reference: 'Functional-unit scaling in life-cycle assessment (ISO 14040 terminology)' },
  { id: 's28-carbon-payback', suite: 28, title: 'Carbon payback time of the hand-worked inventory', mode: 'lca',
    inputs: { ...LCA_HAND, eSave: 50 }, kpi: 'payback', expected: GWP_HAND / (50 * 0.35), tol: 1e-3, kind: 'analytical',
    reference: 'Payback = production footprint/(annual electricity saved × grid carbon intensity)' },
  { id: 's28-mci-linear', suite: 28, title: 'Material circularity indicator of a fully linear product', mode: 'lca',
    inputs: { mode: 'lca', recIn: 0, recOut: 0, recEff: 80, lifeRatio: 1 }, kpi: 'mci', expected: 0.1, tol: 1e-9, kind: 'analytical',
    reference: 'Ellen MacArthur Foundation, Circularity Indicators methodology — MCI = 1 − LFI·0.9/X with LFI = 1' },
  { id: 's28-mci-circular', suite: 28, title: 'Material circularity indicator of a fully circular product', mode: 'lca',
    inputs: { mode: 'lca', recIn: 100, recOut: 100, recEff: 100, lifeRatio: 1 }, kpi: 'mci', expected: 1, tol: 1e-9, kind: 'analytical',
    reference: 'Ellen MacArthur Foundation, Circularity Indicators methodology — LFI = 0 for fully recycled feedstock and lossless recycling' },
  { id: 's28-mci-hand-worked', suite: 28, title: 'Material circularity indicator, hand-worked case (50 % recycled in, 50 % collected, 80 % efficiency)', mode: 'lca',
    inputs: { mode: 'lca', recIn: 50, recOut: 50, recEff: 80, lifeRatio: 1 }, kpi: 'mci', expected: MCI.MCI, tol: 1e-9, kind: 'analytical',
    reference: 'Ellen MacArthur Foundation, Circularity Indicators methodology — V = 0.5, W_C = 0.1, W_F = 0.125, W = 0.6125, LFI = 1.1125/2.0125',
    note: 'LFI = 0.5528, MCI = 0.5025; recycled feedstock is taken at the same process efficiency as end-of-life recycling.' },
  { id: 's28-lfi-hand-worked', suite: 28, title: 'Linear flow index of the hand-worked case', mode: 'lca',
    inputs: { mode: 'lca', recIn: 50, recOut: 50, recEff: 80, lifeRatio: 1 }, kpi: 'lfi', expected: 1.1125 / 2.0125, tol: 1e-9, kind: 'analytical',
    reference: 'Ellen MacArthur Foundation, Circularity Indicators methodology — LFI = (V + W)/(2M + (W_F − W_C)/2)' },
  { id: 's28-mci-lifetime', suite: 28, title: 'Material circularity indicator of a linear product with twice the average lifetime', mode: 'lca',
    inputs: { mode: 'lca', recIn: 0, recOut: 0, recEff: 80, lifeRatio: 2 }, kpi: 'mci', expected: 0.55, tol: 1e-9, kind: 'analytical',
    reference: 'Ellen MacArthur Foundation, Circularity Indicators methodology — utility function F(X) = 0.9/X' },

  // ------------------------------------------------------------------ learning curve, scale economics, make-or-buy
  { id: 's28-wright-exponent', suite: 28, title: 'Wright’s-law exponent for a 15 % learning rate', mode: 'learning',
    inputs: { mode: 'learning', lr: 15 }, kpi: 'b', expected: -0.2345, abs: 5e-5, kind: 'analytical',
    reference: 'Wright’s law, b = log₂(1 − LR) = log₂ 0.85', note: 'Literal −0.2345 (four decimals).' },
  { id: 's28-wright-two-doublings', suite: 28, title: 'Experience curve: cost after two doublings at 20 % learning', mode: 'learning',
    inputs: { mode: 'learning', lr: 20, c0: 100, q0: 1000, qT: 4000 }, kpi: 'c_target_q', expected: 64, tol: 1e-9, kind: 'analytical',
    reference: 'Wright’s law: each doubling of cumulative output multiplies cost by 1 − LR' },
  { id: 's28-wright-target-volume', suite: 28, title: 'Experience curve: cumulative volume to reach 0.8³ of the present cost', mode: 'learning',
    inputs: { mode: 'learning', lr: 20, c0: 100, q0: 1000, cTarget: 51.2 }, kpi: 'q_star', expected: 8000, tol: 1e-9, kind: 'analytical',
    reference: 'Wright’s law inverted, Q* = Q₀(C_target/C₀)^(1/b)', note: 'Three doublings.' },
  { id: 's28-wright-time-to-target', suite: 28, title: 'Years to reach the target with output doubling every year', mode: 'learning',
    inputs: { mode: 'learning', lr: 20, c0: 100, q0: 1000, cTarget: 51.2, qYear: 1000, growth: 100 }, kpi: 't_star', expected: 3, tol: 1e-9, kind: 'analytical',
    reference: 'Geometric series of annual output: Q(t) = Q₀ + q[(1 + g)ᵗ − 1]/g', note: '1000 + 2000 + 4000 kg in three years supplies the 7000 kg needed.' },
  { id: 's28-wright-average-cost', suite: 28, title: 'Average cost over the next production at 50 % learning (logarithmic case)', mode: 'learning',
    inputs: { mode: 'learning', lr: 50, c0: 100, q0: 1000, qT: 4000 }, kpi: 'avg', expected: 100 * 1000 * Math.log(4) / 3000, tol: 1e-9, kind: 'analytical',
    reference: 'Integral of C₀(Q/Q₀)⁻¹ between Q₀ and Q_T divided by the volume produced' },
  { id: 's28-scale-capital-per-unit', suite: 28, title: 'Six-tenths rule: capital cost per unit at tenfold scale', mode: 'learning',
    inputs: { mode: 'learning', sRatio: 10, scaleExp: 0.6 }, kpi: 'cap_ratio', expected: Math.pow(10, -0.4), tol: 1e-9, kind: 'analytical',
    reference: 'Cost–capacity relation: specific capital cost ∝ S^(n−1) (Peters and Timmerhaus)', note: '10^(−0.4) = 0.398.' },
  { id: 's28-make-or-buy', suite: 28, title: 'Make-or-buy break-even volume', mode: 'learning',
    inputs: { mode: 'learning', fixedMake: 250000, varMake: 20, buyPrice: 55 }, kpi: 'q_be', expected: 250000 / 35, tol: 1e-9, kind: 'analytical',
    reference: 'Linear break-even analysis, Q_BE = F/(p_buy − v)' },

  // ------------------------------------------------------------------ material passport and end of life
  { id: 's28-incineration-co2', suite: 28, title: 'CO₂ from complete combustion of an acrylate photopolymer (C₅H₈O₂)', mode: 'passport',
    inputs: { mode: 'passport', route: 'resin', lhv: 0 }, kpi: 'e_inc', expected: 5 * M_CO2 / M_MMA, tol: 2e-3, kind: 'analytical',
    reference: 'Combustion stoichiometry C₅H₈O₂ + 6 O₂ → 5 CO₂ + 4 H₂O with IUPAC standard atomic weights', note: '2.198 kg CO₂/kg; the suite uses the rounded mass ratio 44/12 for CO₂/C.' },
  { id: 's28-recycling-net-cost', suite: 28, title: 'Net recycling cost of ZnO with the recovered zinc valued at the metal price', mode: 'passport',
    inputs: { mode: 'passport', route: 'zno', costRec: 2, valRec: 0, recEff: 80, zinc: 2800 }, kpi: 'c_rec', expected: 2 - 0.8 * (AW.Zn / M_ZNO) * 2.8, tol: 2e-3, kind: 'analytical',
    reference: 'Zinc mass fraction of ZnO (80.34 %) from standard atomic weights; net cost = process cost − efficiency × recovered value' },
  { id: 's28-recycling-net-emissions', suite: 28, title: 'Net emissions of recycling with avoided primary production', mode: 'passport',
    inputs: { mode: 'passport', route: 'zno', eRec: 3, grid: 350, recEff: 80, co2Primary: 12 }, kpi: 'e_rec', expected: 3 * 0.35 - 0.8 * 12, tol: 1e-9, kind: 'analytical',
    reference: 'Avoided-burden (substitution) approach to recycling credits in life-cycle assessment', note: '1.05 − 9.6 = −8.55 kg CO₂e/kg.' },
  { id: 's28-recovery-rate', suite: 28, title: 'Overall material recovery rate', mode: 'passport',
    inputs: { mode: 'passport', recOut: 30, recEff: 80 }, kpi: 'rec_rate', expected: 0.24, tol: 1e-9, kind: 'conservation',
    reference: 'Recovery rate = collection rate × recycling process efficiency' },
];
