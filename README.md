# Designer Materials

Physics-based design, simulation, fabrication, optimisation and reverse engineering of **nanomaterials** and **metamaterials**, organised as 29 connected engineering suites. It runs entirely in the browser, installs on computers, Android phones and iPhones, and keeps working offline.

## What it does

- **29 engineering suites** — from requirements, composition and crystal structure through thermodynamics, nucleation and growth, transport, thermal, electrical, electromagnetic, magnetic, mechanical and acoustic physics, to multiscale homogenisation, control, fabrication, lifetime, metrology, verification/validation, inverse design, scale-up, safety, economics and reverse engineering.
- **One project data bus** — every suite publishes its results; downstream suites pick them up, so a laboratory or industrial case can be run end-to-end from the *Case pipeline* page.
- **Explained outputs** — every number, chart, map and table carries an explanation; every run lists the governing equations it solved and, where chemistry is involved, the balanced chemical, complete-ionic and net-ionic equations and half-cell reactions.
- **Recommendations** — engineering, laboratory, industrial, economic and sustainability actions derived from each run, alongside current literature fetched live.
- **Sensitivity and uncertainty** on any suite: tornado, Morris screening, Sobol indices and Monte Carlo propagation.
- **Materials library** — nanomaterials and metamaterials, including negative-index, negative-permittivity, negative-permeability, auxetic, negative-thermal-expansion, negative-compressibility, negative-stiffness, negative-differential-resistance, negative-effective-mass and negative-differential-thermal-resistance classes.
- **146 creation environments** with the chemical, thermal, transport, field and mechanical conditions each requires.
- **Geometry Studio** — four geometry engines (atomic/crystal, CAD/surface, implicit/lattice, mesh/voxel/image) behind one interface: import from established formats, generate parametrically, or obtain geometry by inverse design. Simulated growth returns geometry that evolves with simulation time.
- **Data portal** — CSV, TSV, instrument exports, Excel, JSON, XML, YAML, Touchstone and more.
- **Credibility gates** — results are labelled by the evidence behind them. A converged calculation is never presented as experimentally validated.

## Honest scope

The solvers are analytical, reduced-order and small-grid numerical models chosen to run on any device, including phones. They are engineering predictions for design, screening and planning; they do not replace calibrated high-fidelity simulation or experiment. Where a high-fidelity external solver is the right tool, the platform exports ready-to-run input decks and geometry. Proprietary CAD formats are recognised and inventoried, with guidance on the neutral export to use.

## Run locally

No build step and no dependencies.

```sh
python3 -m http.server 8080      # then open http://localhost:8080
node tests/run.mjs               # contract, robustness and kernel tests
node tools/build-assets.mjs      # refresh the offline asset list after changing files
```

## Live data

Reference data is requested by the user's browser directly from public services (OpenAlex, Crossref, Europe PMC, Zenodo, DataCite, PubChem, OPTIMADE providers, NOMAD, IMF commodity prices via DBnomics, ECB rates, grid carbon intensity, World Bank) and cached on the device. No server belonging to this project is involved, so freshness never depends on any one machine.

## Hosting and mirrors

The app is a folder of static files. Publish the repository root on any static host; list each address in `mirrors.json` so users can switch if one host is unavailable. Installed copies continue to work offline regardless of host availability. `tools/deploy-mirrors.sh` pushes the same commit to every configured git remote.

## Layout

```
index.html, sw.js, manifest.webmanifest   application shell and offline support
js/core/        numerical kernel, chemistry, plotting, 3-D viewer, sensitivity, live data, engine
js/suites/      s01.js … s29.js — one module per engineering suite
js/geometry/    atomic, implicit/lattice, mesh/CAD, layout engines, format readers and writers
js/pages/       user interface
data/           specification knowledge base, materials library, environments
tests/          automated tests
tools/          build scripts
```
