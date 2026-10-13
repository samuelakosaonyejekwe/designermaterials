# Changes

## 1.4.0
- Published experiments extended to 102 measured results, now covering auxetic lattices, negative effective mass, negative permittivity and permeability, negative-expansion composites, lattice moduli, vacancy concentrations, reliability data, NIST reference materials and certified data-reduction sets. Experiments the solvers do not reproduce are listed in the benchmark notes rather than omitted.
- CATIA V5 parts: the display mesh stored in the file is read (a visualisation mesh at the saved tolerance, not the exact surfaces).
- Autodesk Inventor and Fusion assemblies are placed using the transforms stored in the files; Inventor 2026 files are read.
- Clearer import reports: each file states in one line what was recovered and the best next step.
- The Internet Archive mirror is the single-file edition, served unmodified.

## 1.3.1
- Single-file edition: the complete application as one self-contained document (Install & settings › Mirrors). It runs from any folder, drive or server with no other files, and is the copy kept on the Internet Archive, so that mirror is complete and consistent after a single save.
- Full-resolution runs fall back to computing on the page where background workers are not allowed.

## 1.3.0
- Published experiments: solvers are now compared with measured results from the literature (sources cited with DOI), and each suite's Credibility tab shows whether your case lies near a published experiment — independent evidence without measurements of your own.
- Models corrected where experiments disagreed (results change):
  - Suite 7: ZnO nanorod growth law redesigned around supersaturation at the buffered bath pH with competing bulk precipitation; calibrated on one paper and checked against three others (7–26 % low on the independent set; lengths level off near 1 µm as measured).
  - Suite 12, 5: planar printed split-ring resonator model (4.66 GHz against 4.845 GHz measured for the classic geometry; the previous cylinder formula was ten times too high for printed rings); wire-medium formula valid for thicker wires.
  - Suite 10: spectral phonon model for nanowires and thin films (70 nm silicon wire 27.6 against 28 W/m/K measured; previously 93).
  - Suite 11: electrolyte conductivity with concentration correction (KCl within 1 %).
  - Suite 18: particulate composites at high filler loading (differential scheme as default).
  - Suite 3: zero-point correction for rare-gas solids. Suite 23: BET fitting window chosen by the Rouquerol criteria.
- Geometry: ACIS, Autodesk and JT readers; all configurations of a SOLIDWORKS part; batch-export script generator for closed CAD systems.
- The Internet Archive copy is brought up to date incrementally and verified by content digest on every scheduled run; the schedule keeps itself active.

## 1.2.1
- Parasolid (.x_t, .x_b), SOLIDWORKS part and Siemens NX part geometry is now read: the Parasolid B-rep is translated to STEP and tessellated by the OpenCASCADE kernel.
- A second live address on independent infrastructure (Internet Archive), refreshed automatically at every release together with the content-delivery copy and the source archive.
- The portable launcher fetches the newest server-refreshed data snapshot.

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
