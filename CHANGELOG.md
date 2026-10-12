# Changes

## 1.2.0
- Reference data is now also collected on a server schedule and delivered to every device as one snapshot: live values, the newest arXiv preprints per suite, and measured crystal structures from the Crystallography Open Database.
- Validate and calibrate a case directly on each suite's Credibility tab: enter your measurements to obtain an "experimentally validated" status for that operating point, or fit an input so the model reproduces a measured result.
- Validation benchmarks run in parallel in the background and their results are remembered, so the page opens instantly after the first run.
- The CAD geometry kernels are stored automatically for offline use after installation (skipped when the device asks to save data).
- The portable launcher now starts the newest published release automatically, with the release it was saved from as fallback.
- Benchmark tolerances tightened by improving the models: silicon carrier statistics (temperature-dependent band gap), silver optical constants (measured table), porous-media surface area and permeability estimators, ensemble-averaged statistical models. Related default results in Suites 4, 5, 9, 11, 12, 16, 17 and 25 changed slightly.
- The La-doped SrTiO₃ preset now uses a published Schottky formation energy (1.64 eV per defect; Liu et al., Phys. Chem. Chem. Phys. 16, 15590, 2014).

## 1.1.x
- Independent reference benchmarks for all 29 suites, with a Validation benchmarks page.
- Full-resolution background solver for the hydrothermal ZnO reactor (case-study mesh 100 × 100 × 10, Δt 0.1 s).
- STEP, IGES and BREP import through the OpenCASCADE kernel; Rhino 3DM, IFC and FreeCAD import.
- Portable launcher and source archive.

### Results that changed after verification against reference values
Corrections made while benchmarking alter some results compared with 1.0. If you recorded numbers from 1.0, re-run those cases.
- Suite 3: specific surface area of sheets, slabs and nanotubes (graphene now 2628 m²/g; previously inflated by the vacuum gap); equation-of-state bulk modulus (fit window narrowed).
- Suite 5 and 25: relative density of strut lattices with node overlap; sign of the imaginary effective mass of the mass-in-mass cell.
- Suite 12: negative-index S-parameter retrieval (branch tracking); aluminium Drude constants.
- Suite 14: snap-through extrema of the bistable beam (exact instead of grid-sampled).
- Suite 16: pinned grain growth can no longer exceed the Zener limit.
- Suite 17: Lennard-Jones argon well depth (ε/k_B = 119.8 K).
- Suite 18: bulk modulus of the Kelvin cell.
- Suite 19: finite re-entrant and non-square lattice panels (default 8 × 5 auxetic panel: E 108.3 → 95.5 MPa, ν −0.70 → −0.61).
- Suite 20: end-of-batch horizon of the growth-tracking predictive controller.
- Suite 21: CVD depletion integral near the leading edge.
- Suite 22: sign of the thermal-mismatch stress in the delamination model (default total stress 250 → 350 MPa).
- Suites 23, 24, 26: Student-t quantiles at 1–2 degrees of freedom.
- Suite 29: starting value of the resonator fit.
