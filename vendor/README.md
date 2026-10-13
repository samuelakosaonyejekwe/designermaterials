# Vendored geometry kernels

Third-party WebAssembly geometry kernels used by `js/geometry/kernels.js`. They are **on-demand packs**: nothing in this
directory is loaded when the application starts. A pack is fetched (and then kept by the offline cache) the first time a
file of a matching format is imported. Everything was taken unmodified from the npm registry with `npm pack`, except for
the two loader patches listed below. No `node_modules`, no package dependencies.

| Pack | Package and version | Engine | Formats | Licence | Runtime files | Bytes (runtime) | Bytes (with licence texts) |
|---|---|---|---|---|---|---|---|
| `occt-import-js/` | [occt-import-js](https://www.npmjs.com/package/occt-import-js) 0.0.23 — <https://github.com/kovacsv/occt-import-js> | Open CASCADE Technology — <https://github.com/Open-Cascade-SAS/OCCT> | STEP, IGES, BREP (and the BREP shapes inside FreeCAD documents) | LGPL-2.1 (occt-import-js); LGPL-2.1 with the Open CASCADE exception 1.0 (OCCT) | `occt-import-js.js`, `occt-import-js.wasm` | 7 699 271 | 7 753 900 |
| `rhino3dm/` | [rhino3dm](https://www.npmjs.com/package/rhino3dm) 8.35.0 — <https://github.com/mcneel/rhino3dm> | openNURBS — <https://github.com/mcneel/opennurbs> | Rhino 3DM | MIT (rhino3dm); openNURBS licence (permissive, includes zlib) | `rhino3dm.module.js`, `rhino3dm.wasm` | 2 795 907 | 2 799 181 |
| `web-ifc/` | [web-ifc](https://www.npmjs.com/package/web-ifc) 0.0.78 — <https://github.com/ThatOpen/engine_web-ifc> | web-ifc | IFC | MPL-2.0 | `web-ifc-api.js`, `web-ifc.wasm` | 7 504 318 | 7 521 043 |

Total: 17 999 496 bytes of runtime files (18 074 124 bytes with the licence texts).

Only the single-threaded browser builds are vendored. Not copied: TypeScript declarations, the Node/CommonJS and IIFE
builds, the multi-threaded web-ifc binary (it needs cross-origin isolation), the 64-bit rhino3dm build, worker stubs,
examples and tests.

## Redistribution

Each licence permits redistribution in a public repository, subject to the conditions noted here. This is a summary
of the licence texts kept beside the files, not legal advice.

- **occt-import-js / Open CASCADE Technology (LGPL-2.1).** Redistribution in object form is permitted provided the
  licence text accompanies it (`LICENSE.occt-import-js.txt`, `LICENSE.occt.txt`, `OCCT_LGPL_EXCEPTION.txt`), the
  library can be replaced by the recipient, and the corresponding source is available: it is the upstream source at
  the tagged versions linked above (occt-import-js 0.0.23 and the OCCT revision it pins as a submodule). The pack is a
  separate file that the application loads at run time through its public interface, so it can be replaced by another
  build without touching the application. Supporting documentation must state that the application makes use of Open
  CASCADE Technology — this file and the import report do.
- **rhino3dm (MIT) / openNURBS.** Redistribution permitted with the copyright and permission notice (`LICENSE`,
  `LICENSE.opennurbs.txt`).
- **web-ifc (MPL-2.0).** Redistribution permitted; the file-level copyleft applies to the web-ifc files only. The
  source form of the vendored files is the upstream repository at v0.0.78; `LICENSE.md` is included.

## Modifications

The WebAssembly binaries are byte-identical to the npm distribution. Two JavaScript loaders were changed, and carry a
notice at the top of the file saying so (changed on 2026-10-12):

1. `occt-import-js/occt-import-js.js` and `rhino3dm/rhino3dm.module.js` — Emscripten's embind glue built its call
   wrappers from strings (`new Function(...)`), which a Content-Security-Policy without `'unsafe-eval'` forbids. The
   two factories involved (`craftInvokerFunction` and the `emval` method caller) were replaced by the equivalent
   closure-based implementations that Emscripten itself emits when built with `-sDYNAMIC_EXECUTION=0`; the unused
   string-building helpers (`createJsInvoker`, `newFunc`) were removed. Behaviour is otherwise unchanged.
2. `occt-import-js/occt-import-js.js` only — the package ships a classic script. `export default occtimportjs;` was
   appended so that it can be loaded with a dynamic `import()`, and the two CommonJS `require` calls of the Node branch
   now take the `fs` and `path` modules from the options object (`nodeFs`, `nodePath`), because `require` does not
   exist in an ES module. In a browser that branch is never reached.

`web-ifc/web-ifc-api.js` is unmodified (it is already an ES module and contains no dynamic code generation).

## Content-Security-Policy

The kernels need `script-src 'self' 'wasm-unsafe-eval'` — `'wasm-unsafe-eval'` allows WebAssembly compilation only; no
`'unsafe-eval'` is required. `connect-src 'self'` covers fetching the `.wasm` files. No file here contacts any other
host.

## Updating

Run `npm pack <package>@<version>` in a scratch directory, copy the files named in the table, re-apply the changes
described under "Modifications" (the anchors are the function names given there; verify afterwards that neither
`new Function` nor `eval(` remains), update the versions and byte counts in this file and in the `KERNELS` table of
`js/geometry/kernels.js`, and run `node tests/run.mjs` — `tests/kernels.test.mjs` checks the byte counts, the absence of
dynamic code generation and the analytic volumes.

## fzstd

- Package: `fzstd` 0.1.1 (pure JavaScript Zstandard decompressor), unmodified ES-module build `esm/index.mjs`.
- Licence: MIT (see `fzstd/LICENSE`).
- Upstream: https://github.com/101arrowz/fzstd
- Used for: Zstandard-compressed entries in recent Autodesk Fusion archives.
