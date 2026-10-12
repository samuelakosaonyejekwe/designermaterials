// Geometry import/export hub: format registry, signature-based recognition, parsing, validation, conservative
// repair, suitability assessment, representation conversion and export. Everything is local; procedural files are
// parsed, never executed. Exact B-rep formats (STEP, IGES, BREP, FreeCAD), Rhino 3DM and IFC are tessellated by the
// on-demand WebAssembly kernels in kernels.js, with the built-in lightweight readers as the fallback. Parasolid B-reps
// (transmit files and the partitions inside SOLIDWORKS, NX and Solid Edge documents) are translated to STEP by xt.js and
// tessellated by the same OpenCASCADE kernel.
import * as A from './atomic.js';
import * as F from './formats.js';
import { tpmsVoxels, buildLattice, porousMedium, marchingCubes, latticeToMesh, voxelStats } from './implicit.js';
import { meshStats, repairMesh, bboxOf, triangulate, extrude, mergeMeshes, rodsToMesh, nanorodArray, reactor } from './mesh.js';
import { layoutStats, polyArea, resonator, array as layoutArray } from './layout.js';
import * as K from './kernels.js';
import * as N from './native.js';
import * as XT from './xt.js';

const { min, max, abs, floor, round, hypot, PI } = Math;
const latin = (b, o, n) => { let s = ''; for (let i = o, e = min(b.length, o + n); i < e; i++) s += String.fromCharCode(b[i]); return s; };

// ---------------------------------------------------------------- registry
const SUITES = { cad: [4, 5, 14, 19, 21], mesh: [5, 9, 10, 14, 15, 19], atomic: [3, 17], crystal: [3, 17, 23], layout: [5, 12, 21], voxel: [9, 10, 16, 18, 23, 29], graph: [5, 14, 18], procedural: [5, 25], pointcloud: [23, 29], toolpath: [21, 26], trajectory: [7, 17], orientation: [16, 23], stack: [12, 21], scene: [19] };
const X = {
  cad: 'Export from the authoring CAD system as STEP AP242 with tessellation (or AP214 faceted B-rep), or as STL / 3MF for a surface mesh.',
  mesh: 'Convert with the solver pre-processor (or a mesh converter) to Gmsh MSH 2.2 ASCII, VTU or Abaqus INP.',
  scene: 'Export as glTF / GLB, OBJ, PLY or STL.',
  vox: 'Export as a TIFF stack, NRRD, MetaImage (.mha), VTI or an HDF5 file with a contiguous or gzip-chunked dataset.',
  atom: 'Convert (for example with ASE or Open Babel) to extended XYZ, CIF, POSCAR or PDB.',
  lay: 'Export from the layout tool as GDSII or OASIS (or DXF).',
  pc: 'Export as PLY, PTS, XYZ text or ASCII PCD.',
  traj: 'Convert to an extended-XYZ or LAMMPS text trajectory.',
  proc: 'Run the definition in its own application and export the result: STEP / STL for solids, MSH for meshes, GDSII for layouts, JSON node–edge data for lattices.',
};
// id | name | extensions and special file names | doc sections | representation | priority | support | note | aliases (;) | advice key
const TABLE = `
step|STEP / ISO 10303|.step .stp .p21 .stpz|1,7,8,9,10,11,13,15,16,17|cad|P1|full|With the on-demand OpenCASCADE kernel (vendor/occt-import-js, about 7.7 MB, fetched on first use): every B-rep face (planar, analytic, swept, trimmed NURBS) of every solid, shell and assembly instance is tessellated to a triangle mesh in millimetres with body names and colours, at a linear deflection of 0.001 × the mean bounding-box edge and an angular deflection of 0.5 rad unless set otherwise. The result is a tessellation of the exact B-rep: the surface definitions themselves, PMI and the assembly tree are not retained. Gzip- or zip-wrapped files (.stpZ, .stp.gz) are unpacked first. Without the kernel (offline before it has been cached) the built-in reader gives header, products, units and entity inventory, meshes faceted B-rep, AP242 tessellation and planar straight-edged faces, and shows curved faces as a wireframe only.|STEP;STEP assemblies;NURBS geometry|cad
iges|IGES|.iges .igs|1,7,8,10,11,13|cad|P1|full|With the on-demand OpenCASCADE kernel (vendor/occt-import-js, about 7.7 MB, fetched on first use): trimmed and untrimmed surfaces (types 108, 114, 118, 120, 122, 128, 143, 144 …) and manifold solids (186) are tessellated to a triangle mesh in millimetres at the stated deflection; surface models whose faces do not share edges come through as open (unsewn) shells and are reported as such. The exact surfaces are not retained. Without the kernel the built-in reader gives the global section (units, author), the entity inventory and a wireframe of points, lines, arcs and polylines.||cad
brep|OpenCASCADE BREP|.brep .brp|1,2|cad|P2|full|OpenCASCADE's native shape text format, tessellated by the on-demand OpenCASCADE kernel (vendor/occt-import-js, about 7.7 MB, fetched on first use) exactly as STEP is (all faces, at the stated deflection; the exact surfaces are not retained). A BREP file declares no length unit and carries no names or colours. Without the kernel: section inventory only.|BREP;BRP|cad
parasolid|Parasolid|.x_t .x_b .xmt_txt .xmt_bin|1,8|cad|P1|partial|Text and neutral-binary transmit files (parts, assemblies of instanced bodies, partitions) are read node by node through the schema the file embeds, translated entity by entity to STEP AP214 and tessellated by the on-demand OpenCASCADE kernel (vendor/occt-import-js, about 7.7 MB, fetched on first use) in millimetres. Written exactly: planes, cylinders, cones, spheres, tori, NURBS, swept, spun and offset surfaces; lines, circles, ellipses and NURBS curves. Approximated and counted in the report: edges on intersection curves and surface (SP) curves (cubic B-splines through the stored chart or sampled points) and rolling-ball blends between analytic surfaces (B-spline surfaces through circular sections). Not translated (the face is left out, the body becomes an open shell without a volume): cliff-edge blends, blends on NURBS or swept surfaces, foreign (PE) geometry, facet (mesh) bodies; wire bodies and mirrored instances are skipped. Measured on 66 public transmit files and partitions (Parasolid 19 to 36, written by SOLIDWORKS, NX and other Parasolid applications): 8 899 of 8 963 faces translated (99.3 %). Names of bodies are kept; feature history, sketches, PMI and assembly constraints are not in a transmit file. Schemas older than version 12 and bare-binary files are not read. Without the kernel: header only.||cad
sat|ACIS SAT|.sat|1|cad|P1|partial|Header (version, product, date, units scale), entity inventory, and the topological vertex points with edges drawn as straight chords (a wireframe preview with a vertex bounding box, scaled to millimetres). Surfaces and curves are not evaluated: no faces are tessellated and body transforms are not applied.|SAT|cad
sab|ACIS SAB|.sab|1|cad|P2|inventory|Recognised by the binary ACIS signature; not decoded.||cad
jt|JT|.jt|1|cad|P2|inventory|Header and table of contents: version, byte order and the count and size of every segment type (scene graph, shape LODs, JT / XT B-rep, PMI, meta data). The entropy-coded tessellation and B-rep segments are not decoded.||cad
dxf|DXF|.dxf|1,2,8,9,10,11,15|layout|P1|partial|ASCII DXF: LINE, LWPOLYLINE/POLYLINE (bulges), CIRCLE, ARC, ELLIPSE, SOLID, 3DFACE; blocks, splines, hatches and text are counted, not converted.|AutoCAD|cad
dwg|DWG|.dwg|1,2|cad|P3|inventory|Release identified from the AC10xx signature; the proprietary object database is not decoded.|BricsCAD|Save as DXF (ASCII) for 2-D geometry, or export STEP for solids.
ifc|IFC|.ifc|1|cad|P3|partial|With the on-demand web-ifc kernel (vendor/web-ifc, about 7.5 MB): the body geometry of every building element (swept and extruded solids, B-reps, Boolean clippings and openings, mapped items) is tessellated in metres, one body per element with its IFC class, name and colour; curved profiles use a fixed number of segments per circle. Property sets, materials, the spatial structure and relationships are not imported, and representation types the kernel does not evaluate are simply absent. Without the kernel: header and counts of building elements only.||Export the equipment or structure of interest as STEP, or as a tessellated glTF / OBJ / STL.
3dm|Rhino 3DM|.3dm|1,2,13|cad|P2|partial|With the on-demand openNURBS kernel (vendor/rhino3dm, about 2.8 MB): mesh objects in full; Breps and extrusions through the render meshes Rhino saved in the file; untrimmed NURBS surfaces by sampling; SubD as its control net; block instances placed (one level); units, layers, names and colours kept. openNURBS has no surface mesher, so trimmed Breps saved without render meshes are counted but not tessellated. Without the kernel: archive version only.|Rhino;3DM|Save from Rhino with render meshes (shaded view, not "Save small"), or export STEP, to bring trimmed NURBS solids through.
dae|COLLADA|.dae|1|scene|P3|inventory|XML scene recognised; geometry libraries counted.||scene
gltf|glTF|.gltf .glb|1|scene|P2|partial|Triangle primitives with node transforms from GLB or embedded buffers; external buffers, materials, animation not imported.||scene
usd|USD|.usd .usda .usdc .usdz|1|scene|P2|inventory|Text / crate / package variant identified; prim types counted for .usda.||scene
vrml|VRML|.wrl .vrml|1|scene|P3|partial|IndexedFaceSet geometry only; Transform nodes are not applied.||scene
x3d|X3D|.x3d .x3db|1|scene|P3|partial|XML IndexedFaceSet geometry only; binary .x3db and Transform nodes not handled.||scene
solidworks|SOLIDWORKS|.sldprt .sldasm .slddrw|2|cad|P2|partial|Both the OLE2 layout (2014 and earlier) and the 2015+ compressed-stream container are opened: preview image, document and custom properties, release, configuration names, feature names and the stream inventory are recovered, and the Parasolid partition of each part configuration is unpacked. The partition of the first configuration is translated to STEP and tessellated by the on-demand OpenCASCADE kernel (vendor/occt-import-js, about 7.7 MB, fetched on first use) as described under Parasolid. Measured on the eleven NIST CTC / FTC parts (SOLIDWORKS 2018) against the published STEP files: 2 312 of 2 330 faces translated, nine of eleven bodies closed, and volume within 0.17 %, area within 0.01 % and bounding box within 0.005 % for every part whose published STEP is itself a closed solid; the two open parts lose faces on blends against NURBS surfaces. Tolerant models can leave sliver gaps that keep a body from closing. The feature tree, sketches, dimensions, PMI, mates and display lists are not interpreted; assemblies and drawings hold no partition and give the inventory only. Without the kernel: inventory and the unpacked partition.|SLDPRT;SLDASM|cad
catia|CATIA V5|.catpart .catproduct .catdrawing|2|cad|P2|inventory|V5 container signature recognised; not decoded.|CATPart|cad
nx|Siemens NX|.prt|2|cad|P2|partial|The OLE2 layout (NX 12 and earlier) and the SPLMSSTR container (NX 1847 and later) are recognised and the Parasolid partition is located by signature (stored raw or deflated), translated to STEP and tessellated by the on-demand OpenCASCADE kernel (vendor/occt-import-js, about 7.7 MB, fetched on first use) as described under Parasolid. Measured on the eleven NIST CTC / FTC parts (NX 1980) and six NIST D2MI parts (NX 9): 2 735 of 2 737 faces translated and ten of the eleven CTC / FTC bodies closed; against the published STEP files the volume is within 0.12 % (most within 0.005 %) and the area within 0.2 % wherever the published STEP is a closed solid, and the three D2MI parts compared agree to better than 0.002 %. Sketch, section and reference bodies stored as separate transmit streams and roll-back deltas are left aside. Features, expressions, sketches, PMI, drawings and assembly structure are not interpreted (an assembly file holds no component geometry). A .prt file without either container signature is only assumed to be NX. Without the kernel: container inventory and the located partition.|NX PRT|cad
creo|PTC Creo|.prt .asm|2|cad|P3|inventory|Recognised by the #UGC header; release line read. Not decoded.||cad
inventor|Autodesk Inventor|.ipt .iam|2|cad|P2|partial|OLE2 container opened: thumbnail, summary and design-tracking properties (part number, material, designer, dates, saving release) and the segment inventory are recovered. The model segments (ShapeManager B-rep, graphics) are in Inventor's proprietary encoding and are not decoded: no geometry.|Inventor IPT|cad
fusion|Autodesk Fusion|.f3d .f3z|2|cad|P3|inventory|Zip archive: entry listing, preview image when one is stored, nested component archives and the location and size of the ShapeManager B-rep segments (.smb / .smbh). Those segments are an ACIS derivative that no open kernel reads: no geometry.||cad
solidedge|Solid Edge|.par .psm .asm|2|cad|P3|inventory|OLE2 container opened: property sets, thumbnail when present and the stream inventory. Streams are searched for Parasolid transmit data by signature and, when any is found, it is translated and tessellated as described under Parasolid; this has not been verified against real Solid Edge files, whose model streams may use an encoding the search does not recognise.||cad
fcstd|FreeCAD|.fcstd|2|cad|P3|partial|Zip archive: document properties, object list and thumbnail from Document.xml; with the on-demand OpenCASCADE kernel (vendor/occt-import-js, about 7.7 MB, fetched on first use) the OpenCASCADE BREP shape of every visible object is tessellated in millimetres (hidden Boolean inputs and features already contained in a Body are left out, App::Part placements applied, object labels and colours kept). App::Link and assembly instances are not expanded and the parametric history is not imported. Without the kernel: object inventory only.||cad
onshape|Onshape||2|cad|P3|inventory|Cloud-native: there is no file to open.||Use Onshape's export to STEP or Parasolid text, then STEP here.
revit|Autodesk Revit|.rvt .rfa|2|cad|P3|inventory|OLE2 container opened: property sets, preview image when present and the stream inventory. The model streams are not decoded.||Export IFC or, for single equipment items, STEP / STL.
stl|Stereolithography mesh|.stl|3,4,7,8,9,10,11,15,17|mesh|P1|full|ASCII and binary (decided by the size rule); facets welded to an indexed mesh. No units, colours or materials.|STL;STL-based surface mesh|
obj|Wavefront OBJ|.obj|3,7,8,15|mesh|P1|full|Vertices, faces (polygons fan-triangulated), line elements and group names; no materials or textures.|OBJ|
ply|Polygon File Format|.ply|3,7,14,15|mesh|P1|full|ASCII and binary (little/big endian); faces or bare point clouds; colours and normals dropped.|PLY;Stanford triangle format|
3mf|3D Manufacturing Format|.3mf|3,8,11,15|mesh|P1|full|Core-specification meshes with build-item and component transforms and units; materials, beam-lattice and slice extensions are reported as not imported.|3MF|
amf|Additive Manufacturing File Format|.amf|3,8,15|mesh|P1|partial|Plain or zipped XML: vertices and triangles of all objects with units; curved triangles, materials and constellations not applied.|AMF|
off|OFF|.off|3|mesh|P1|full|OFF / COFF / NOFF with polygon fan-triangulation; colours dropped.|COFF|
gts|GTS|.gts|3|mesh|P3|full|GNU Triangulated Surface text format.||
3ds|3DS|.3ds|3|scene|P3|inventory|Chunk signature recognised; not decoded.||scene
fbx|FBX|.fbx|3|scene|P3|inventory|Binary or ASCII variant and version identified; not decoded.||scene
abc|Alembic|.abc|3,19|scene|P2|inventory|Ogawa or HDF5 container recognised; time-sampled geometry not decoded.||Export the frames of interest as OBJ / PLY, or as a .pvd + .vtu / .vtp series.
vdb|OpenVDB|.vdb|3,7,8,11,12,13,19|voxel|P2|inventory|Magic number and file version read; sparse trees are not decoded.|OpenVDB sequences|Resample to a dense grid and export as VTI, NRRD or a TIFF stack (for example from Houdini or Blender).
msh|Gmsh|.msh|4,8,10,11,17|mesh|P1|partial|MSH 2.2 and 4.0/4.1 ASCII: nodes, elements (orders above one reduced to corners) and physical/entity tags; binary files are not read.|MSH;Gmsh MSH;multiregion volume meshes|mesh
vtk|VTK legacy|.vtk|4,7,8,11,12,13,16,17|mesh|P1|partial|Legacy ASCII POLYDATA, UNSTRUCTURED_GRID and STRUCTURED_POINTS; legacy BINARY and structured/rectilinear grids are not read.|VTK|mesh
vtu|VTK Unstructured Grid|.vtu|4,10,16|mesh|P1|full|ASCII, inline base64 and appended data, optionally zlib-compressed (LZ4 / LZMA are not supported); field arrays are listed, not attached.|VTU|
vtp|VTK PolyData|.vtp|4|mesh|P1|full|Same encodings as VTU; polygons, strips and lines.||
vtm|VTK Multiblock|.vtm|4|mesh|P2|inventory|Block list with referenced files; open the referenced .vtu / .vtp / .vti individually.||Load the referenced block files directly.
ensight|EnSight|.case .geo|4|mesh|P3|inventory|Case file: geometry and variable file names; the geometry file itself is not decoded.||mesh
cgns|CGNS|.cgns|4|mesh|P2|inventory|HDF5 container: zone and array inventory (names, shapes); zone connectivity is not assembled.||mesh
unv|Universal File|.unv|4|mesh|P1|full|Datasets 2411 (nodes), 2412 (elements), 164 (units) and groups.|UNV|
bdf|NASTRAN Bulk Data|.bdf .nas .dat|4,11,18|mesh|P1|partial|GRID, CTRIA3/6, CQUAD4/8, CTETRA, CHEXA, CPENTA, CPYRAM, CBAR/CBEAM/CROD in small, large and free field; local coordinate systems, properties and loads not applied.|NASTRAN BDF;BDF;NASTRAN-style meshes where supported|mesh
inp|Abaqus Input|.inp|4,11,18|mesh|P1|partial|*NODE, *ELEMENT (continuum, shell, beam/truss families), *NSET/*ELSET sizes; instance placements, materials and steps not applied.|Abaqus INP;INP|mesh
cdb|ANSYS Mechanical APDL|.cdb|4|mesh|P2|inventory|NBLOCK / EBLOCK sizes and element-type declarations; blocks are not decoded.|CDB|mesh
lsdyna|LS-DYNA Keyword|.k .key|4|mesh|P2|inventory|Keyword inventory (node and element card counts, parts); not decoded.|LS-DYNA KEY|mesh
exodus|Exodus II|.e .exo .ex2|4,19|mesh|P2|inventory|netCDF classic or HDF5 container recognised; variables listed for the HDF5 flavour.||mesh
xdmf|XDMF|.xdmf .xmf|4,17,19|mesh|P2|partial|Grid, time-step and heavy-data inventory; inline XML grids, and HDF5-backed grids when the .h5 file is supplied, are reconstructed (first grid only).|XDMF + HDF5|mesh
med|MED|.med|4|mesh|P2|inventory|HDF5 container: mesh and field inventory; not assembled.||mesh
su2|SU2|.su2|4|mesh|P2|full|Points, elements and boundary markers of the native ASCII mesh.||
fluent|Fluent mesh|.msh .msh.h5|4|mesh|P3|inventory|ASCII case/mesh: dimension, node/face/cell counts and zone names; face-based connectivity is not assembled. HDF5 variant: inventory.||mesh
openfoam|OpenFOAM|polyMesh points faces owner neighbour boundary|4|mesh|P2|inventory|A polyMesh is a directory of several files; a single file gives only its FoamFile header and entry count.|OpenFOAM polyMesh|Run foamToVTK and import the .vtu / .vtk, or convert with the solver tools to MSH.
plot3d|PLOT3D|.xyz .p3d .q|4|mesh|P3|inventory|Structured multi-block grid: block dimensions from an ASCII header when present; coordinates not decoded.||mesh
xyz|XYZ atomic coordinates|.xyz|5,6,7,19|atomic|P1|full|Element symbols and Cartesian coordinates in ångström; multi-frame files give animation frames.|XYZ|
extxyz|Extended XYZ|.extxyz .xyz|5,6,7,19|atomic|P1|full|Lattice, per-frame comment properties and trajectories.|EXTXYZ;Extended XYZ trajectory|
pdb|Protein Data Bank|.pdb|5,7|atomic|P1|full|ATOM / HETATM records and CRYST1 cell; connectivity, anisotropic factors and biological assemblies not applied.|PDB|
mmcif|PDBx/mmCIF|.cif .mmcif|5|atomic|P2|inventory|Recognised by the _atom_site.Cartn categories; entry id and atom count reported.||Convert to PDB or extended XYZ (for example with gemmi, ASE or Open Babel).
mol|MDL MOL|.mol|5|atomic|P1|full|V2000 atom and bond blocks.|MOL|
sdf|Structure Data File|.sdf|5|atomic|P1|full|V2000 records; the first molecule is the geometry, further records are counted.|SDF|
mol2|MOL2|.mol2|5|atomic|P2|full|Tripos atoms, bonds and optional CRYSIN cell.|Tripos MOL2|
cml|Chemical Markup Language|.cml|5|atomic|P3|full|atomArray / bondArray with x3, y3, z3 coordinates.||
gaussian|Gaussian input|.gjf .com|5|atomic|P3|partial|Cartesian molecule specification; Z-matrix input is not converted.||atom
zmat|Z-matrix|.zmat .zmt|5|atomic|P3|inventory|Internal coordinates recognised; not converted to Cartesian.||atom
psf|Protein Structure File|.psf|5|atomic|P3|inventory|Topology only (no coordinates): atom, bond, angle counts.||Load the matching coordinate file (PDB / CRD) instead.
charmmcrd|CHARMM coordinates|.crd .cor|5|atomic|P3|inventory|Atom count and title lines; element assignment from CHARMM atom types is not attempted.||atom
gro|GROMACS structure|.gro|5|atomic|P2|full|Coordinates (nm → Å) and box vectors.|GRO|
gmxtop|GROMACS topology|.top .itp|5|atomic|P3|inventory|Force-field topology (no coordinates): molecule types and counts.||Load the matching .gro / .pdb coordinate file.
lammpsdata|LAMMPS data|.data .lmp|5,6,7|atomic|P2|full|atomic / charge / full styles with orthogonal or triclinic box; type → element mapping uses Masses when present.|LAMMPS data/dump;LAMMPS atomic structure|
lammpsdump|LAMMPS dump|.lammpstrj .dump|5,6,19|trajectory|P2|full|Text dumps: first frame as geometry, further frames (capped) as animation.|LAMMPS trajectories|
pqr|PQR|.pqr|5|atomic|P3|partial|Read through the PDB reader (whitespace-delimited charge / radius columns are ignored).||atom
poscar|VASP POSCAR|POSCAR CONTCAR .vasp .poscar|5,6,7,16|crystal|P1|full|VASP 4/5, direct or Cartesian, selective dynamics.|VASP POSCAR / CONTCAR;POSCAR;CONTCAR|
xsf|XSF|.xsf .axsf|5,6|crystal|P1|full|PRIMVEC / PRIMCOORD / ATOMS blocks; grids and forces dropped.|XCrySDen Structure Format|
cif|Crystallographic Information File|.cif|5,6,7,10,16|crystal|P1|full|Cell, atom sites and symmetry operators when listed; occupancy and displacement parameters are kept only where the atomic engine stores them.|CIF|
xdatcar|VASP XDATCAR|XDATCAR|6,19|trajectory|P2|partial|Lattice and configurations with a constant cell: first frame as geometry, further frames (capped at 40) as animation.||
qein|Quantum ESPRESSO input|.in .pwi|6|crystal|P2|partial|CELL_PARAMETERS / ATOMIC_POSITIONS of pw.x input; ibrav-defined lattices other than those the atomic engine supports are rejected with a message.|Quantum ESPRESSO|atom
qeout|Quantum ESPRESSO output|.out .pwo|6|crystal|P2|inventory|Program banner, atom count and final energy line; relaxed coordinates are not extracted.||atom
wien2k|WIEN2k structure|.struct|6|crystal|P3|inventory|Title and lattice line only.||atom
fhiaims|FHI-aims geometry|geometry.in|6|crystal|P3|partial|atom / atom_frac / lattice_vector lines.||atom
abinit|ABINIT input|.abi .in|6|crystal|P3|inventory|Variable inventory (natom, ntypat, acell); not converted.||atom
cp2k|CP2K input|.inp|6|crystal|P2|partial|&CELL (A, B, C or ABC) and &COORD blocks; SCALED coordinates and external coordinate files are not followed.|CP2K geometry inputs|atom
castep|CASTEP cell|.cell|6|crystal|P2|full|LATTICE_CART / LATTICE_ABC and POSITIONS_FRAC / POSITIONS_ABS blocks.|CASTEP CELL|
shelx|SHELX|.ins .res|6|crystal|P3|inventory|TITL, CELL and SFAC lines reported; symmetry expansion is not performed.||Convert to CIF (for example with shelxl, Olex2 or cif tools).
jana|JANA structure data|.m40 .m50|6|crystal|P3|inventory|Recognised by extension; modulated-structure parameters are not interpreted.|specialised crystallographic structure files|Export the average structure as CIF.
crystalmaker|CrystalMaker document|.cmdf|6|crystal|P3|inventory|Proprietary document; recognised by extension only.||Export as CIF or XYZ.
gds|GDSII / GDS|.gds .gdsii .gds2|7,8,9,10,16|layout|P1|full|BOUNDARY, BOX, PATH (outlined), SREF / AREF hierarchy (flattened), units, layers / datatypes and text labels.|GDSII|
oasis|OASIS|.oas .oasis|8,9,10,16|layout|P1|full|RECTANGLE, POLYGON, PATH, TRAPEZOID, CTRAPEZOID, CIRCLE, PLACEMENT with all repetition types, and DEFLATE CBLOCKs; properties and text are skipped.||
svg|SVG|.svg|8,9,11,15|layout|P2|partial|path (lines, Béziers, arcs flattened), polygon, polyline, rect, circle, ellipse, line with group transforms; no fill rules, clipping, <use> or text.||lay
caltechcif|CIF / Caltech Intermediate Form|.cif|9|layout|P3|inventory|Recognised by DS / L / B / P commands; counts of boxes, polygons, wires, layers and symbol definitions.||lay
gerber|Gerber|.gbr .ger .gtl .gbl|9,10|layout|P3|partial|RS-274X: standard apertures (C, R, O, P), linear and circular draws, flashes and regions; aperture macros, clear-polarity subtraction and step-repeat are not applied.||lay
excellon|Excellon drill|.drl .xln|9|layout|P3|partial|Tool table and drill hits as circles; routed slots as single hits.|Excellon|
lef|LEF|.lef|9|layout|P3|inventory|Units, layer and macro names.||lay
def|DEF|.def|9|layout|P3|inventory|Design name, units, die area, component and net counts.||lay
openaccess|OpenAccess||9|layout|P3|inventory|A database, not a file: nothing to parse.||Stream the cell view out as GDSII or OASIS from the EDA tool.
odbpp|ODB++||9|layout|P3|inventory|Directory / archive of many files; only the archive listing can be shown.||Export Gerber RS-274X + Excellon, or DXF.
tiff|TIFF / BigTIFF|.tif .tiff|7,11,12,17,20|voxel|P1|full|Classic and BigTIFF, strips or tiles, none / LZW / PackBits / Deflate, 1/8/16/32-bit and float, grey or RGB; pages become z slices. JPEG-compressed and planar-RGB files are rejected with a message.|TIFF stacks;EBSD image maps;voxel TIFF;voxel volumes;voxel data;voxel models;voxel formats|
ometiff|OME-TIFF|.ome.tif .ome.tiff|12|voxel|P2|partial|Pixel data as for TIFF with voxel size from the OME-XML header; channel / time dimensions are not separated.||vox
hdf5|HDF5|.h5 .hdf5 .hdf .he5|7,8,10,12,13,14,16,17,18,19|voxel|P1|partial|Classic layout: old-style and compact-link groups, contiguous / compact / chunked (deflate, shuffle) numeric datasets; the largest 3-D dataset becomes the volume. Dense groups, new chunk indices and other filters are inventoried only.|HDF5 graph datasets;Signed Distance Field;Level-set geometry|vox
nifti|NIfTI|.nii .nii.gz|12|voxel|P2|full|NIfTI-1 single file (optionally gzip); intensity scaling applied; orientation matrices not applied.||
dicom|DICOM|.dcm|12|voxel|P3|inventory|Part-10 preamble recognised and a few header tags read; pixel data (often compressed, one slice per file) are not decoded.||Convert the series to NIfTI, NRRD or a TIFF stack (for example with dcm2niix or Fiji).
nrrd|NRRD|.nrrd .nhdr|12|voxel|P2|full|Raw, gzip and ASCII encodings; detached headers need the data file supplied alongside.||
metaimage|MetaImage|.mhd .mha|12|voxel|P3|full|Local or companion raw data, optional zlib compression.||
vti|VTK ImageData|.vti|12,13|voxel|P1|full|Cell or point scalars in ASCII, base64 or appended form (zlib optional).|VTI|
raw|RAW binary volume|.raw .bin .vol|12|voxel|P3|partial|Headerless: dimensions, sample type and byte order must be supplied (a perfect-cube byte count is the only guess).||vox
mrc|MRC / CCP4|.mrc .map .ccp4 .rec|12|voxel|P2|full|Modes 0, 1, 2, 6 with cell-derived voxel size; non-standard axis order is reported, not permuted.|MRC|
em|EM format|.em|12|voxel|P3|partial|512-byte header (dimensions, sample type) followed by raw data.||vox
zarr|Zarr|.zarr .zarray|12|voxel|P3|inventory|A chunked directory store: array metadata (.zarray) can be shown, chunks are separate files.||vox
omezarr|OME-Zarr||12|voxel|P2|inventory|Multiscale Zarr hierarchy: directory store, not importable as one file.||vox
scad|OpenSCAD|.scad|13|procedural|P1|partial|Parsed, not executed: cube, sphere, cylinder, translate / rotate / scale / mirror, union / difference / intersection with numeric variables; evaluated as a signed-distance field and meshed. Modules, loops, extrusions, hull and minkowski are reported as unsupported.|Constructive Solid Geometry|proc
python|Python geometry scripts|.py|13|procedural|P2|inventory|Never executed: imports, functions and geometry-library usage are listed.|Parametric CAD scripts;Python parameterisations|proc
grasshopper|Grasshopper definitions|.gh .ghx|13|procedural|P2|inventory|.ghx XML: component count and names; binary .gh recognised by extension. The definition cannot be evaluated without Rhino.|GH/GHX|proc
dynamo|Dynamo definitions|.dyn|13|procedural|P3|inventory|JSON / XML graph: node count and types; not evaluated.||proc
geo|Gmsh geometry scripts|.geo|13|procedural|P1|partial|Parsed, not executed: variables, Point, Line, Circle, Spline as a wireframe; loops, surfaces, volumes, physical groups and OpenCASCADE primitives are inventoried.|Gmsh GEO|proc
json|JSON|.json|10,13,16,18|procedural|P1|partial|This application's geometry export, plain geometry objects, node–edge lattices, vertex–face meshes, generator definitions (TPMS, lattice, porous medium, resonator array, nanorod array, reactor) and layer-stack tables. Arbitrary schemas are not guessed.|JSON graph definitions;JSON-based unit cells;Lattice structure data;Custom unit-cell descriptors;Implicit mathematical geometry;implicit equations;procedural geometry;procedural lattice definitions;implicit mathematical definitions;procedural definitions;implicit geometry;assemblies and custom parameterised definitions;parametric definitions;Layer-stack tables;procedural descriptions;statistical RVE data|
yaml|YAML|.yaml .yml|13,16|procedural|P3|partial|Simple block YAML mapped onto the same schemas as JSON.||
xml|XML|.xml|13,16|stack|P3|inventory|Generic XML: root element and tag counts; <layer material thickness> tables are read as a layer stack.||Use a format with a defined geometry schema (JSON node–edge data, 3MF, VTU, GDSII).
csv|CSV|.csv .tsv|14,16,17,18|pointcloud|P2|partial|Tables with x, y, z columns (point cloud), node / edge or segment tables (lattice), or material / thickness columns (layer stack). Crease patterns and other graphs are read as plain node–edge networks (fold assignments are not interpreted).|CSV node/edge tables;node–edge CSV;CSV coordinates;Node–edge data;beam/node data;Graph data;graph-based representations;crease-pattern data;crease graphs|
mtex|MTEX-compatible orientation tables|.txt|20|orientation|P3|inventory|Column table with Euler-angle headers recognised; columns and row count reported.||Export the map as CTF or ANG to keep the scan grid.
xyzpoints|XYZ point cloud|.xyz .asc|14|pointcloud|P3|full|Whitespace- or comma-separated x y z rows (extra columns dropped).||
pts|PTS|.pts|14|pointcloud|P3|full|Count line followed by x y z [intensity r g b] rows.||
ptx|PTX|.ptx|14|pointcloud|P3|partial|First scan block; the registration matrix is not applied.||pc
las|LAS|.las|14|pointcloud|P3|inventory|Public header: version, point count, scale, offset and bounding box; point records are not decoded.||pc
laz|LAZ|.laz|14|pointcloud|P3|inventory|As LAS (header only); points are LASzip-compressed.||pc
e57|E57|.e57|14|pointcloud|P3|inventory|ASTM E57 signature and version; XML section and compressed vectors not decoded.||pc
pcd|PCD|.pcd|14|pointcloud|P3|partial|ASCII and uncompressed binary; binary_compressed is rejected.||pc
cli|CLI|.cli|15|toolpath|P2|partial|ASCII Common Layer Interface: layers and closed contours; hatches summarised; binary CLI rejected.||Export ASCII CLI or the source STL / 3MF.
slc|SLC|.slc|15|toolpath|P2|inventory|ASCII header keywords (units, extents, type); binary contour data not decoded.||Export ASCII CLI or the source STL / 3MF.
gcode|G-code|.gcode .nc .gco .ngc|15|toolpath|P2|partial|G0–G3 motion with units and absolute/relative modes: work moves as polylines plus layer, length and time statistics; machine-specific codes ignored.||
buildfile|Build-file formats||15|toolpath|P3|inventory|Machine-vendor job files are proprietary.|Proprietary project and manufacturing databases|Export the neutral input of the build (3MF / STL geometry, CLI slices or G-code).
graphml|GraphML|.graphml|18|graph|P2|full|Nodes with x / y / z data keys and edges (optional radius attribute); nodes without coordinates are placed on a circle and flagged.||
gexf|GEXF|.gexf|18|graph|P2|full|Nodes with viz:position and edges.||
gml|GML|.gml|18|graph|P3|full|Graph Modelling Language nodes (x, y, z or graphics block) and edges.||
pvd|VTK time series|.pvd|19|trajectory|P2|partial|Time-step list; referenced .vtu / .vtp / .vti frames are loaded when supplied alongside (capped at 40).|PVD|
xtc|GROMACS trajectories|.xtc .trr|19|trajectory|P2|inventory|First-frame header (atom count, step, time); compressed coordinates not decoded and no topology is present.|XTC;TRR|traj
ctf|Oxford EBSD CTF|.ctf|20|orientation|P2|full|Header, phases and the phase / Euler-angle / band-contrast map on the scan grid.|CTF;EBSD maps|
ang|EDAX / TSL EBSD|.ang|20|orientation|P2|full|Header and the phase / Euler-angle / image-quality map; hexagonal grids are resampled to a rectangular array.|ANG|
h5ebsd|HDF5-based EBSD|.h5ebsd .h5oina .oh5|20|orientation|P2|inventory|HDF5 inventory (scan groups, array shapes); vendor layouts are not mapped to a grid.|H5EBSD|Export the scan as CTF or ANG.
dream3d|DREAM.3D data|.dream3d|17,20|orientation|P2|partial|HDF5 inventory; a Phases (or small-label) cell array is imported as a voxel volume when readable. FeatureIds, orientations and statistics are not imported.|DREAM3D;Dream3D-style datasets|vox
`;
export const FORMATS = TABLE.trim().split('\n').map(l => {
  const [id, name, ex, fam, representation, priority, support, note, al, adv] = l.split('|'), t = ex.split(/\s+/).filter(Boolean);
  return { id, name, ext: t.filter(e => e[0] === '.'), names: t.filter(e => e[0] !== '.'), family: representation, families: fam.split(',').map(Number), representation, priority, support, suites: SUITES[representation] || [], note, aliases: al ? al.split(';') : [], neutral: adv ? (X[adv] || adv) : '' };
});
const UNKNOWN = { id: 'unknown', name: 'Unrecognised file', ext: [], names: [], family: 'unknown', families: [], representation: 'scene', priority: 'P3', support: 'inventory', suites: [], note: 'Neither the signature, the structure nor the extension matches a known geometry format.', aliases: [], neutral: 'Export to one of the supported neutral formats: STEP, STL, 3MF, MSH, VTU, CIF, extended XYZ, GDSII, TIFF stack or HDF5.' };
const BY = Object.fromEntries(FORMATS.map(f => [f.id, f]));
const EXT = new Map(); for (const f of FORMATS) for (const e of f.ext) { if (!EXT.has(e)) EXT.set(e, []); EXT.get(e).push(f.id); }
export const UNKNOWN_FORMAT = UNKNOWN;
const normName = s => String(s).toLowerCase().replace(/[^a-z0-9+]/g, '');
// Registry search by format name, alias, id or extension (case and punctuation insensitive).
export function lookupFormats(term) { const k = normName(term); return k ? FORMATS.filter(f => normName(f.name) === k || f.id === k || f.aliases.some(a => normName(a) === k) || f.ext.some(e => e.slice(1) === k) || f.names.some(n => normName(n) === k)) : []; }
// detectFormat may also return the pseudo-ids 'unknown', 'gzip', 'compress' (.Z) and 'zipwrap' (a zip holding one CAD file): wrappers that importFile unpacks.
export const findFormat = id => BY[id] || (['unknown', 'gzip', 'compress', 'zipwrap'].includes(id) ? UNKNOWN : null);

// ---------------------------------------------------------------- detection
const MAGIC = [ // [offset, bytes-as-latin1, id, description]
  [0, '%SEMI-OASIS\r\n', 'oasis', 'OASIS magic string'], [0, 'glTF', 'gltf', 'binary glTF magic'], [0, 'AutoCAD Binary DXF', 'dxf', 'binary DXF sentinel'], [0, 'V5_CFV2', 'catia', 'CATIA V5 container signature'],
  [0, '#UGC:', 'creo', 'Creo "#UGC" header'], [0, 'ACIS BinaryFile', 'sab', 'ACIS binary signature'], [0, '3D Geometry File Format', '3dm', 'Rhino openNURBS header'], [0, 'Kaydara FBX Binary', 'fbx', 'binary FBX signature'], [0, '; FBX', 'fbx', 'ASCII FBX header'],
  [0, 'Ogawa', 'abc', 'Alembic Ogawa signature'], [0, ' BDV', 'vdb', 'OpenVDB magic number'], [0, 'PXR-USDC', 'usd', 'USD crate signature'], [0, '#usda', 'usd', 'USD text header'], [128, 'DICM', 'dicom', 'DICOM part-10 prefix'],
  [0, 'LASF', 'las', 'LAS public header signature'], [0, 'ASTM-E57', 'e57', 'E57 signature'], [208, 'MAP ', 'mrc', '"MAP " stamp at byte 208'], [0, 'NRRD000', 'nrrd', 'NRRD magic'], [0, '-SLCVER', 'slc', 'SLC header keyword'], [0, 'C Binary', 'ensight', 'EnSight binary geometry header'],
  [0, 'PS\0\0\0', 'parasolid', 'Parasolid binary transmit header'], [0, 'SPLMSSTR', 'nx', 'Siemens NX part container signature (SPLMSSTR)'], [0, 'Channel Text File', 'ctf', 'Oxford CTF first line'], [0, 'PSF', 'psf', 'PSF keyword'],
];
const TEXT_SIG = [ // strong structural signatures on the first 8 kB of text
  ['step', /^\s*ISO-10303-21\s*;/], ['dxf', /^\s*0\s*\r?\n\s*SECTION\s*\r?\n|^\s*999\s*\r?\n/], ['ply', /^ply\r?\n/], ['vtk', /^# vtk DataFile Version/], ['msh', /^\s*\$MeshFormat/], ['xdmf', /<Xdmf\b/], ['x3d', /<X3D\b/], ['dae', /<COLLADA\b/],
  ['amf', /<amf\b/], ['graphml', /<graphml\b/i], ['gexf', /<gexf\b/i], ['svg', /<svg\b/], ['cml', /<(cml|molecule)\b[^>]*>/], ['grasshopper', /<Archive name="Root">/], ['unv', /^\s{4}-1\s*\r?\n\s*\d+\s*\r?\n/], ['lsdyna', /^\*KEYWORD\b/m],
  ['cdb', /^\/(COM|PREP7|NOPR)\b|^NBLOCK,|^EBLOCK,/m], ['ensight', /^FORMAT\s*\r?\n\s*type:\s*ensight/im], ['openfoam', /FoamFile\s*\{/], ['su2', /^\s*NDIME\s*=/m], ['fluent', /^\s*\(\s*[012]\s+["(\d]/], ['cli', /\$\$HEADERSTART/],
  ['gerber', /%FS[LT][AI]X\d\dY\d\d\*%|%MO(MM|IN)\*%|%ADD\d+[A-Za-z]/], ['excellon', /^M48\s*$/m], ['mhd', /^\s*(ObjectType\s*=\s*Image|NDims\s*=)/m], ['pcd', /^(#\s*\.PCD|VERSION\s+[\d.]+\s*\r?\nFIELDS)/], ['ang', /^#\s*(TEM_PIXperUM|x-star|WorkingDistance|GRID:)/m],
  ['mol2', /@<TRIPOS>MOLECULE/], ['lammpsdump', /^ITEM:\s*TIMESTEP/], ['xsf', /^\s*(CRYSTAL|SLAB|POLYMER|MOLECULE|ANIMSTEPS\s+\d+)\s*\r?\n[\s\S]*?(PRIMVEC|ATOMS|PRIMCOORD)|^\s*ATOMS\s*\r?\n/], ['qeout', /Program PWSCF/], ['castep', /%BLOCK\s+(LATTICE_|POSITIONS_)/i],
  ['cp2k', /^\s*&(GLOBAL|FORCE_EVAL|MOTION)\b/im], ['gmxtop', /^\s*\[\s*(moleculetype|defaults|atomtypes)\s*\]/m], ['jt', /^Version \d+\.\d+ JT/], ['brep', /^\s*(DBRep_DrawableShape\s+)?CASCADE Topology V\d/], ['parasolid', /^\*\*ABCDEFGHIJKLMNOPQRSTUVWXYZ[\s\S]{0,400}\*\*PARASOLID|^\*\*PART1;/], ['lammpsdata', /^\s*\d+\s+atoms\s*$[\s\S]*\bxlo\s+xhi\b/m],
];
const isNum = s => s !== '' && Number.isFinite(Number(s));
const ELEM = /^(?:[A-Z][a-z]?\d*[+-]?|\d{1,3})$/;
function xyzKind(head) { // atomic XYZ / extended XYZ / point cloud / PLOT3D
  const L = head.split(/\r?\n/), l0 = (L[0] || '').trim(), t2 = (L[2] || '').trim().split(/\s+/);
  if (/^\d+$/.test(l0) && t2.length >= 4 && ELEM.test(t2[0]) && isNum(t2[1]) && isNum(t2[2]) && isNum(t2[3])) return /Lattice\s*=|Properties\s*=/i.test(L[1] || '') ? ['extxyz', 'atom count, then a comment line with Lattice= / Properties=, then element rows'] : ['xyz', 'atom-count line followed by element + x y z rows'];
  const rows = L.filter(l => l.trim() && !/^[#%]/.test(l.trim())).slice(0, 12).map(l => l.trim().split(/[\s,;]+/));
  if (rows.length > 1 && rows[0].length <= 3 && rows[0].every(s => /^\d+$/.test(s)) && rows.slice(1).every(r => r.every(isNum)) && rows.slice(1, 4).some(r => r.length !== 3 && r.length !== 4 && r.length !== 6 && r.length !== 7)) return ['plot3d', 'integer block-dimension header followed by unstructured runs of numbers'];
  if (rows.length && rows.slice(rows[0].every(isNum) ? 0 : 1).every(r => r.length >= 3 && isNum(r[0]) && isNum(r[1]) && isNum(r[2]))) return ['xyzpoints', 'purely numeric x y z rows without element symbols'];
  return null;
}
function cifKind(t) {
  if (/_atom_site\.(Cartn_x|group_PDB|label_atom_id)|^_entry\.id\b|_pdbx_/m.test(t)) return ['mmcif', 'data_ block with dotted _atom_site.* / _entry.* categories (PDBx/mmCIF)'];
  if (/^\s*data_/m.test(t) && /_cell_length_a|_atom_site_fract_x|_atom_site_label|_symmetry_|_space_group_/.test(t)) return ['cif', 'data_ block with _cell_length / _atom_site items (crystallographic CIF)'];
  if (/^\s*(DS\s+\d+|L\s+\w+\s*;|B\s+-?\d+\s+-?\d+\s+-?\d+)/m.test(t) && /;\s*(\r?\n|$)/.test(t) && !/^\s*data_/m.test(t)) return ['caltechcif', 'DS / L / B commands terminated by semicolons (Caltech Intermediate Form layout)'];
  if (/^\s*data_/m.test(t)) return ['cif', 'data_ block (CIF syntax) without recognised cell items'];
  return null;
}
const hit = (id, confidence, why) => ({ id, confidence, why });
function extsOf(base) { const lower = base.toLowerCase(), parts = lower.split('.'), out = []; for (let i = 1; i < parts.length; i++) out.push('.' + parts.slice(i).join('.')); return out; }
const WRAPPED = /\.(stp|step|p21|igs|iges|brep|brp|ifc|3dm|stl|obj|ply|x_t|sat|jt)$/i;
function zipKind(b, exts) {
  let names = []; try { names = F.unzip(b, /^$/).map(e => e.name); } catch (e) { return hit('unknown', 0.3, 'zip signature but the central directory could not be read'); }
  const has = re => names.some(n => re.test(n));
  if (has(/3dmodel\.model$/i) || has(/\.model$/i)) return hit('3mf', 0.99, 'zip package containing a 3D model part (OPC / 3MF)');
  if (has(/\.amf$/i) || exts.includes('.amf')) return hit('amf', 0.9, 'zip archive holding an AMF document');
  if (has(/^Document\.xml$/)) return hit('fcstd', 0.97, 'zip archive with Document.xml (FreeCAD)');
  if (has(/\.usd[ac]?$/i) || exts.includes('.usdz')) return hit('usd', 0.9, 'zip package of USD layers (usdz)');
  if (has(/(^|\/)\.zarray$|(^|\/)\.zgroup$|(^|\/)zarr\.json$/)) return hit(has(/\.zattrs$/) ? 'omezarr' : 'zarr', 0.85, 'zipped Zarr store (.zarray / .zgroup entries)');
  if (has(/(^|\/)matrix\/matrix$/i)) return hit('odbpp', 0.85, 'archive with matrix/matrix (ODB++ job)');
  for (const [e, id] of [['.f3d', 'fusion'], ['.f3z', 'fusion'], ['.sldprt', 'solidworks'], ['.sldasm', 'solidworks'], ['.3mf', '3mf']]) if (exts.includes(e)) return hit(id, 0.7, `zip container with extension ${e}`);
  const wrapped = names.filter(n => WRAPPED.test(n)); if (wrapped.length === 1 || (wrapped.length && exts.some(e => ['.stpz', '.zip'].includes(e)))) return hit('zipwrap', 0.9, `zip archive wrapping ${wrapped[0]}`);
  return hit('unknown', 0.3, `zip archive with ${names.length} entries (${names.slice(0, 4).join(', ')}${names.length > 4 ? ', …' : ''}) that match no known geometry package`);
}
// Signature and structure first, extension second. Returns {id, confidence, why}.
export function detectFormat(name, bytes) {
  const b = F.toBytes(bytes ?? ''), base = String(name || '').split(/[\\/]/).pop(), lower = base.toLowerCase(), exts = extsOf(base), ext = exts[exts.length - 1] || '', has = e => exts.includes(e);
  const head = latin(b, 0, 8192), dv = b.length >= 8 ? new DataView(b.buffer, b.byteOffset, b.byteLength) : null;
  let nul = 0; for (let i = 0, e = min(b.length, 4096); i < e; i++) if (b[i] === 0) nul++;
  const utf16 = b.length > 1 && ((b[0] === 0xFF && b[1] === 0xFE) || (b[0] === 0xFE && b[1] === 0xFF)), binary = nul > 0 && !utf16, text = binary ? '' : (utf16 ? F.toText(b, 16384) : head);
  const byName = FORMATS.find(f => f.names.some(n => n.toLowerCase() === lower || (n === n.toUpperCase() && n.length > 4 && lower.startsWith(n.toLowerCase()))));
  // --- binary signatures
  if (head.startsWith('**ABCDEFGHIJKLMNOPQRSTUVWXYZ') && /\*\*PARASOLID !/.test(head)) return hit('parasolid', 0.97, 'Parasolid transmit header (character-set preamble and **PARASOLID record)');
  if (b.length >= 8 && b[0] === 0x89 && head.startsWith('HDF\r\n\x1a\n', 1)) {
    const scan = latin(b, 0, min(b.length, 262144));
    for (const [e, id] of [['.h5ebsd', 'h5ebsd'], ['.h5oina', 'h5ebsd'], ['.oh5', 'h5ebsd'], ['.dream3d', 'dream3d'], ['.cgns', 'cgns'], ['.med', 'med'], ['.msh.h5', 'fluent'], ['.e', 'exodus'], ['.exo', 'exodus'], ['.ex2', 'exodus'], ['.abc', 'abc']]) if (has(e)) return hit(id, 0.95, `HDF5 signature with extension ${e}`);
    if (/CGNSLibraryVersion/.test(scan)) return hit('cgns', 0.9, 'HDF5 file containing a CGNSLibraryVersion node');
    if (/DataContainerBundles|DREAM3D|FileVersion.{0,40}DREAM/.test(scan)) return hit('dream3d', 0.85, 'HDF5 file with DREAM.3D container names');
    if (/INFOS_GENERALES|ENS_MAA/.test(scan)) return hit('med', 0.9, 'HDF5 file with MED group names');
    return hit('hdf5', 0.99, 'HDF5 signature \\x89HDF\\r\\n\\x1a\\n');
  }
  if (head.startsWith('CDF\x01') || head.startsWith('CDF\x02')) return hit('exodus', has('.e') || has('.exo') || has('.ex2') ? 0.9 : 0.5, 'netCDF classic signature (Exodus II is netCDF-based; a generic netCDF file has the same signature)');
  if (head.startsWith('PK\x03\x04')) return zipKind(b, exts);
  if (b.length >= 6 && b[0] === 0 && b[1] === 6 && b[2] === 0 && b[3] === 2) return hit('gds', 0.99, 'GDSII HEADER record (00 06 00 02)');
  if (/^(II\*\0|MM\0\*|II\+\0|MM\0\+)/.test(head)) return hit(has('.ome.tif') || has('.ome.tiff') || /<OME\b/.test(latin(b, 0, 65536)) || /<OME\b/.test(latin(b, max(0, b.length - 65536), 65536)) ? 'ometiff' : 'tiff', 0.99, `TIFF byte-order mark and magic number${head[2] === '+' || head[3] === '+' ? ' (BigTIFF)' : ''}`);
  if (/^AC10\d\d/.test(head)) return hit('dwg', 0.99, `DWG release signature ${head.slice(0, 6)}`);
  if (b.length >= 8 && b[0] === 0xD0 && b[1] === 0xCF && b[2] === 0x11 && b[3] === 0xE0) {
    for (const [e, id] of [['.sldprt', 'solidworks'], ['.sldasm', 'solidworks'], ['.slddrw', 'solidworks'], ['.ipt', 'inventor'], ['.iam', 'inventor'], ['.par', 'solidedge'], ['.psm', 'solidedge'], ['.asm', 'solidedge'], ['.rvt', 'revit'], ['.rfa', 'revit'], ['.prt', 'nx']]) if (has(e)) return hit(id, 0.9, `OLE2 compound file with extension ${e}`);
    const u = latin(b, 0, min(b.length, 131072)).replace(/\0/g, '');
    if (/UG_PART/.test(u)) return hit('nx', 0.8, 'OLE2 compound file with an NX UG_PART stream'); if (/RSeDb|RSeStorage/.test(u)) return hit('inventor', 0.8, 'OLE2 compound file with Inventor stream names'); if (/BasicFileInfo/.test(u)) return hit('revit', 0.8, 'OLE2 compound file with a Revit BasicFileInfo stream'); if (/JSitesList|PARASOLID/.test(u)) return hit('solidedge', 0.6, 'OLE2 compound file with Solid Edge-style streams');
    return hit('unknown', 0.3, 'OLE2 compound file (used by several CAD systems and office documents); the extension does not identify the application');
  }
  if (b.length > 64 && b[4] === 0 && b[5] === 0 && b[6] === 0 && b[7] === 4) { for (let i = 8; i < 60; i++) if (b[i] === 0x14 && b[i + 1] === 0 && b[i + 2] === 6 && b[i + 3] === 0 && b[i + 4] === 8 && b[i + 5] === 0) return hit('solidworks', has('.sldprt') || has('.sldasm') || has('.slddrw') ? 0.97 : 0.85, 'SOLIDWORKS 2015+ container (compressed-stream marker 14 00 06 00 08 00 after the file header)'); }
  if (b[0] === 0x1F && b[1] === 0x9D) return hit('compress', 0.9, 'Unix compress (.Z) container');
  for (const [off, sig, id, why] of MAGIC) if (b.length >= off + sig.length && latin(b, off, sig.length) === sig) { if (id === 'las') return hit(has('.laz') ? 'laz' : 'las', 0.98, why); if (id === 'psf' && !/^PSF\b/.test(head)) continue; return hit(id, 0.98, why); }
  if (b[0] === 0x1F && b[1] === 0x8B) { if (has('.nii.gz')) return hit('nifti', 0.9, 'gzip stream with extension .nii.gz'); if (has('.svgz')) return hit('svg', 0.9, 'gzip stream with extension .svgz'); return hit('gzip', 0.9, 'gzip container'); }
  if (dv && b.length >= 348 && (dv.getInt32(0, true) === 348 || dv.getInt32(0, false) === 348) && /^n[+i]1/.test(latin(b, 344, 3))) return hit('nifti', 0.99, 'NIfTI-1 header (sizeof_hdr = 348, magic n+1 / ni1)');
  if (dv && b.length >= 84) { const n = dv.getUint32(80, true); if (84 + 50 * n === b.length && n > 0) return hit('stl', /^\s*solid\b/.test(head) ? 0.97 : 0.99, `binary STL size rule: 84 + 50 × ${n} facets = ${b.length} bytes${/^\s*solid\b/.test(head) ? ' (the header begins with "solid" but the size rule shows it is binary)' : ''}`); }
  if (dv && binary && dv.getInt32(0, false) === 1995 && (has('.xtc') || b.length > 12)) return hit('xtc', has('.xtc') ? 0.95 : 0.6, 'XTC frame magic number 1995');
  if (dv && binary && dv.getInt32(0, false) === 1993 && has('.trr')) return hit('xtc', 0.9, 'TRR frame magic number 1993');
  if (binary && b[0] === 0x4D && b[1] === 0x4D && has('.3ds')) return hit('3ds', 0.9, '3DS main chunk id 0x4D4D');
  if (binary) {
    if (/^ply\r?\n/.test(head)) return hit('ply', 0.99, 'PLY header (binary payload)'); if (/^(#\s*\.PCD|VERSION)/.test(head) && /DATA\s+binary/.test(head)) return hit('pcd', 0.97, 'PCD header with binary data');
    if (/^NRRD000/.test(head)) return hit('nrrd', 0.99, 'NRRD magic'); if (/ElementDataFile\s*=\s*LOCAL/.test(head)) return hit('metaimage', 0.97, 'MetaImage header with local binary data'); if (/<VTKFile\b/.test(head)) { /* handled below as text */ } else if (/^# vtk DataFile/.test(head)) return hit('vtk', 0.97, 'VTK legacy header (BINARY payload)');
    else {
      if (byName) return hit(byName.id, 0.7, `file name ${base}`);
      const c = exts.map(e => EXT.get(e)).find(Boolean);
      if (c) { if (ext === '.prt') return hit('nx', 0.5, 'binary .prt without the Creo "#UGC" header: assumed Siemens NX (not verifiable from the signature)'); if (ext === '.xyz' || ext === '.q' || ext === '.p3d') return hit('plot3d', 0.6, 'binary file with a PLOT3D extension'); if (ext === '.msh') return hit('unknown', 0.3, 'binary .msh that is neither Gmsh (no $MeshFormat) nor an ASCII Fluent mesh'); if (ext === '.cif') return hit('unknown', 0.3, 'binary data in a .cif file: not a crystallographic, macromolecular or Caltech CIF');
        const id = c.find(i => BY[i].support === 'inventory' || ['raw', 'em', 'mrc', 'stl', 'gds'].includes(i)) || c[0]; return hit(id, 0.55, `binary content; identified by extension ${ext} only`); }
      return hit('unknown', 0.2, 'binary file with no known signature or extension');
    }
  }
  // --- text structure
  if (/^\s*ISO-10303-21\s*;/.test(text)) return hit(/FILE_SCHEMA\s*\(\s*\(\s*'IFC/i.test(text) ? 'ifc' : 'step', 0.99, /FILE_SCHEMA\s*\(\s*\(\s*'IFC/i.test(text) ? 'ISO 10303-21 file whose FILE_SCHEMA is an IFC schema' : 'ISO 10303-21 header');
  const l0 = text.split(/\r?\n/, 1)[0] || '';
  if (l0.length >= 73 && l0[72] === 'S' && /^\s*\d+\s*$/.test(l0.slice(73))) return hit('iges', 0.98, 'IGES fixed-format record: section letter S in column 73');
  const vx = /<VTKFile\b[^>]*\btype\s*=\s*"(\w+)"/.exec(text);
  if (vx) { const id = ({ UnstructuredGrid: 'vtu', PolyData: 'vtp', ImageData: 'vti', vtkMultiBlockDataSet: 'vtm', Collection: 'pvd' })[vx[1]]; return id ? hit(id, 0.99, `VTK XML file of type ${vx[1]}`) : hit(has('.vti') ? 'vti' : 'vtu', 0.5, `VTK XML file of type ${vx[1]}, which this importer does not reconstruct`); }
  if (/^\s*solid\b/.test(text) && /facet\s+normal|endsolid|outer\s+loop/.test(text)) return hit('stl', 0.97, 'ASCII STL: "solid" followed by facet records (size rule for binary not satisfied)');
  if (/^(ST)?C?N?4?n?OFF\b/.test(text)) return hit('off', 0.97, 'OFF header keyword');
  if (/^\s*\d{3,5}\s+\d+\s+\d+\s+\d+\s*\r?\n/.test(text) && (has('.sat') || /ACIS|\bbody\b|\blump\b/.test(text))) return hit('sat', 0.9, 'ACIS SAT header line (version, record counts)');
  if (/V[23]000\s*$/m.test(text.slice(0, 600))) return hit(/^\$\$\$\$/m.test(F.toText(b)) || has('.sdf') ? 'sdf' : 'mol', 0.95, 'MDL connection table counts line (V2000 / V3000)');
  for (const [id, re] of TEXT_SIG) if (re.test(text)) { if (id === 'fluent' && !has('.msh') && !has('.cas')) continue; if (id === 'gerber' && /^\s*</.test(text)) continue; return hit(id === 'mhd' ? 'metaimage' : id, 0.95, `${BY[id === 'mhd' ? 'metaimage' : id].name} structure recognised in the file content`); }
  if (/^\s*(BEGIN\s+BULK|CEND\b|SOL\s+\d|GRID\*?\s*[, ]\s*\d|GRID\*?\s{1,8}\d)/mi.test(text) && !has('.inp')) return hit('bdf', 0.9, 'NASTRAN bulk-data cards (GRID / BEGIN BULK)');
  if (/^\*(HEADING|NODE|PART|ELEMENT)\b/mi.test(text)) return hit('inp', 0.95, 'Abaqus keyword lines (*HEADING / *NODE / *ELEMENT)');
  if (/^(HEADER|CRYST1|ATOM  |HETATM|MODEL )/m.test(text) && /^(ATOM  |HETATM)/m.test(text)) return hit(has('.pqr') ? 'pqr' : 'pdb', 0.95, 'PDB-style ATOM / HETATM records');
  if (/&(control|system|electrons)\b/i.test(text) && /ATOMIC_(POSITIONS|SPECIES)/i.test(F.toText(b))) return hit('qein', 0.95, 'Quantum ESPRESSO namelists and ATOMIC_POSITIONS card');
  if (/^\s*(atom|atom_frac|lattice_vector)\s+-?[\d.]/m.test(text)) return hit('fhiaims', 0.9, 'FHI-aims atom / lattice_vector lines');
  if (/^\s*(DESIGN\s+\S+\s*;)/m.test(text) && /\b(DIEAREA|COMPONENTS|UNITS DISTANCE)\b/.test(text)) return hit('def', 0.9, 'DEF DESIGN / DIEAREA statements');
  if (/^\s*(MACRO\s+\S+|UNITS\s*$|LAYER\s+\S+\s*$)/m.test(text) && /\bEND\s+(LIBRARY|UNITS|\S+)\b/.test(text) && (has('.lef') || /\bSITE\b|\bPIN\b/.test(text))) return hit('lef', 0.85, 'LEF MACRO / LAYER / UNITS blocks');
  if (/^TITL\b/m.test(text) && /^CELL\b/m.test(text)) return hit('shelx', 0.9, 'SHELX TITL and CELL instructions');
  if (has('.cif') || has('.mmcif') || /^\s*data_\S*\s*$/m.test(text)) { const k = cifKind(F.toText(b, 200000)); if (k) return hit(k[0], 0.93, k[1]); }
  if (has('.xyz') || has('.extxyz')) { const k = xyzKind(text); if (k) return hit(k[0], 0.9, `.xyz disambiguated by content: ${k[1]}`); }
  if (/\b(Point|Line|Circle|Curve Loop|Line Loop|Plane Surface)\s*\(\s*[\w+*\- ]+\s*\)\s*=\s*\{/.test(text) || /SetFactory\s*\(/.test(text)) return hit('geo', 0.93, 'Gmsh geometry statements (Point(…) = {…};)');
  if (has('.geo') && /^.*\r?\n.*\r?\nnode id (assign|given|off|ignore)/.test(text)) return hit('ensight', 0.9, 'EnSight Gold ASCII geometry header ("node id …" on line 3)');
  const tr = text.trimStart();
  if (tr[0] === '{' || tr[0] === '[') { let j = null; if (b.length < 6e7) { try { j = JSON.parse(F.toText(b)); } catch (e) { j = null; } } if (j) { if (j.asset && (j.meshes || j.scenes || j.nodes)) return hit('gltf', 0.97, 'JSON with glTF "asset" and mesh / scene arrays'); if ((j.Nodes && j.Connectors) || has('.dyn')) return hit('dynamo', 0.9, 'JSON graph with Nodes and Connectors (Dynamo)'); if (j.zarr_format !== undefined) return hit('zarr', 0.95, 'Zarr array metadata (zarr_format key)'); return hit('json', 0.9, 'valid JSON document'); } }
  if (/^\s*(<\?xml|<!DOCTYPE|<[A-Za-z])/.test(text)) { if (/<Workspace\b[^>]*\bVersion/.test(text) && has('.dyn')) return hit('dynamo', 0.9, 'Dynamo XML workspace'); if (/<model\b[^>]*3dmanufacturing/.test(text)) return hit('unknown', 0.4, 'bare 3MF model XML outside its zip package: re-zip as .3mf'); if (!(has('.gml') && /\bgraph\s*\[/.test(text))) return hit('xml', has('.xml') ? 0.8 : 0.5, 'generic XML document (no known geometry root element)'); }
  if (/\bgraph\s*\[/.test(text) && /\bnode\s*\[/.test(text)) return hit('gml', 0.9, 'Graph Modelling Language "graph [ node [ … ] ]" blocks');
  if (/^#\s/m.test(text) && has('.ang')) return hit('ang', 0.8, '.ang file with a "#" header block');
  // --- file name and extension
  if (byName) return hit(byName.id, 0.8, `file name ${base}`);
  const scadish = /\b(cube|sphere|cylinder|translate|rotate|difference|union|intersection|linear_extrude|module)\s*\(/.test(text), pyish = /^\s*(import\s+\w+|from\s+\w+\s+import|def\s+\w+\s*\(|class\s+\w+)/m.test(text);
  if (has('.scad') || (scadish && !pyish && /;\s*$/m.test(text) && !ext)) return hit('scad', has('.scad') ? 0.9 : 0.6, has('.scad') ? 'extension .scad with OpenSCAD syntax' : 'OpenSCAD module calls');
  if (has('.py') || (pyish && !ext)) return hit('python', 0.85, 'Python source (import / def statements)');
  const csvish = () => { const t = F.parseTable(text.split(/\r?\n/).slice(0, 40).join('\n')), h = t.header || []; if (h.some(c => ['source', 'from', 'target', 'x1', 'node1'].includes(c)) || text.split(/\r?\n/).some(l => F.EDGE_HEADER.test(l))) return hit('csv', 0.85, 'table with node / edge or segment columns'); if (h.some(c => /euler|phi1|^phi$/.test(c))) return hit('mtex', 0.8, 'table with Euler-angle columns'); if (h.includes('x') && h.includes('y')) return hit('csv', 0.8, 'table with x, y(, z) columns'); if (h.some(c => /thick/.test(c))) return hit('csv', 0.8, 'table with a thickness column (layer stack)'); if (!t.header && t.rows.length && t.rows.every(r => r.length >= 3 && r.slice(0, 3).every(isNum))) return hit(has('.csv') ? 'csv' : 'xyzpoints', 0.6, 'numeric table with at least three columns (read as x y z)'); return null; };
  if (ext === '.dat' || ext === '.txt' || ext === '.csv' || ext === '.tsv' || ext === '.asc' || !ext) { const k = csvish(); if (k) return k; if (ext === '.txt' || ext === '.dat' || !ext) return hit('unknown', 0.25, `text file${ext ? ' with the ambiguous extension ' + ext : ''} whose content matches no known geometry structure`); }
  if (ext === '.inp') return hit('unknown', 0.3, '.inp file that is neither an Abaqus keyword file nor a CP2K input (the extension is shared by many codes)');
  if (ext === '.in') { if (/^\s*(acell|rprim|xred|xcart|znucl|ntypat|natom)\b/m.test(text)) return hit('abinit', 0.85, 'ABINIT variables (acell / rprim / xred …)'); return hit('unknown', 0.3, '.in file that is not a Quantum ESPRESSO, ABINIT or FHI-aims geometry input'); }
  if (ext === '.out') return hit('unknown', 0.3, '.out file without the Quantum ESPRESSO banner');
  if (ext === '.msh') return hit('unknown', 0.3, '.msh file with neither a Gmsh $MeshFormat section nor Fluent section syntax');
  if (ext === '.prt') return hit('nx', 0.4, '.prt text without the Creo "#UGC" header');
  if (ext === '.geo') return hit('geo', 0.5, 'extension .geo (no Gmsh entity statements found near the start)');
  if (ext === '.gml') return hit('unknown', 0.3, '.gml file that is not a Graph Modelling Language graph (possibly Geography Markup Language)');
  if (ext === '.nc' || ext === '.gcode' || ext === '.gco' || ext === '.ngc' || /^(\s*(N\d+\s*)?[GM]\d+\b[^\n]*\n){3}/m.test(text)) { if (/^\s*(N\d+\s*)?[GM]\d+/m.test(text)) return hit('gcode', exts.some(e => ['.gcode', '.nc', '.gco', '.ngc'].includes(e)) ? 0.9 : 0.6, 'G / M command words at line starts'); }
  if (ext === '.gro' || (/^\s*\d+\s*$/.test(text.split(/\r?\n/)[1] || '') && /^\s*\d+[A-Za-z]/.test(text.split(/\r?\n/)[2] || '') && !ext)) return hit('gro', ext === '.gro' ? 0.85 : 0.6, 'GROMACS fixed-column structure (title, atom count, residue rows)');
  if (/^\s*[-+]?\d*\.?\d+(?:[eE][-+]?\d+)?\s*\r?\n(\s*[-+\d.eE]+\s+[-+\d.eE]+\s+[-+\d.eE]+\s*\r?\n){3}/.test(text.split(/\r?\n/).slice(1).join('\n')) && !ext) return hit('poscar', 0.7, 'POSCAR layout: comment, scale factor and three lattice-vector lines');
  if (/^\s*#[pntPNT]?\s+\S+/m.test(text) && (has('.gjf') || has('.com'))) return hit('gaussian', 0.85, 'Gaussian route section (# …)');
  for (const e of exts) { const c = EXT.get(e); if (c) { const id = e === '.xyz' ? 'xyz' : e === '.cif' ? 'cif' : c[0]; return hit(id, c.length > 1 ? 0.45 : 0.6, c.length > 1 ? `extension ${e} only (shared by ${c.map(i => BY[i].name).join(', ')}); the content gave no decisive structure` : `extension ${e} only`); } }
  return hit('unknown', 0.1, b.length ? 'no signature, structure or extension matched' : 'empty file');
}

// ---------------------------------------------------------------- small local parsers and header sniffers
const txt = c => c._t ?? (c._t = F.toText(c.bytes));
const atomsOut = (s, notes = [], lost = []) => { const g = { kind: 'atoms', atoms: s.atoms, cell: s.cell || null, title: s.title || '', units: 'Å', notes: [...notes, ...(s.warnings || [])], lost, meta: {} }; for (const k of ['formula', 'spacegroup', 'bonds', 'nFrames', 'charges', 'residues', 'symmetry', 'properties']) if (s[k] !== undefined) g.meta[k] = s[k]; if (s.frames && s.frames.length > 1) g.frames = s.frames.slice(0, 40).map(f => ({ t: f.t, atoms: f.atoms, ...(f.cell ? { cell: f.cell } : {}) })); return g; };
function parseAims(t) { const atoms = [], cell = [], frac = []; for (const l of t.split(/\r?\n/)) { const a = l.replace(/#.*/, '').trim().split(/\s+/); if (a[0] === 'lattice_vector') cell.push(a.slice(1, 4).map(Number)); else if (a[0] === 'atom') atoms.push({ el: a[4], x: +a[1], y: +a[2], z: +a[3] }); else if (a[0] === 'atom_frac') frac.push(a); } if (frac.length) { if (cell.length !== 3) throw new Error('FHI-aims: atom_frac lines need three lattice_vector lines'); for (const a of frac) { const u = [+a[1], +a[2], +a[3]]; atoms.push({ el: a[4], x: u[0] * cell[0][0] + u[1] * cell[1][0] + u[2] * cell[2][0], y: u[0] * cell[0][1] + u[1] * cell[1][1] + u[2] * cell[2][1], z: u[0] * cell[0][2] + u[1] * cell[1][2] + u[2] * cell[2][2], f: u }); } } if (!atoms.length || atoms.some(a => !a.el || !Number.isFinite(a.x + a.y + a.z))) throw new Error('FHI-aims: no valid "atom x y z Element" lines'); return { atoms, cell: cell.length === 3 ? cell : null, title: 'FHI-aims geometry' }; }
function parseCP2K(t) {
  const blk = n => { const m = new RegExp('&' + n + '\\b[^\\n]*\\n([\\s\\S]*?)&END\\s+' + n, 'i').exec(t); return m ? m[1] : null; }, co = blk('COORD'), ce = blk('CELL'); if (!co) throw new Error('CP2K: no &COORD block (coordinates may be in an external file named by COORD_FILE_NAME)');
  let cell = null; if (ce) { const v = k => { const m = new RegExp('^\\s*' + k + '\\s+(?:\\[\\w+\\]\\s+)?([-\\d.eE+]+)\\s+([-\\d.eE+]+)\\s+([-\\d.eE+]+)', 'mi').exec(ce); return m ? [+m[1], +m[2], +m[3]] : null; }, a = v('A'), b = v('B'), c = v('C'), abc = v('ABC'); if (a && b && c) cell = [a, b, c]; else if (abc) cell = [[abc[0], 0, 0], [0, abc[1], 0], [0, 0, abc[2]]]; }
  if (/^\s*SCALED\b(?!\s+(F|\.FALSE\.))/mi.test(co)) throw new Error('CP2K: SCALED (fractional) coordinates are not supported by this reader');
  const atoms = []; for (const l of co.split(/\r?\n/)) { const a = l.replace(/[!#].*/, '').trim().split(/\s+/); if (a.length >= 4 && /^[A-Za-z]/.test(a[0]) && a.slice(1, 4).every(isNum)) atoms.push({ el: a[0].replace(/[\d_].*$/, '').replace(/^(.)(.?)/, (m, p, q) => p.toUpperCase() + q.toLowerCase()), x: +a[1], y: +a[2], z: +a[3] }); }
  if (!atoms.length) throw new Error('CP2K: &COORD block holds no "Element x y z" rows'); return { atoms, cell, title: 'CP2K input' };
}
function parseXDATCAR(t) {
  const L = t.split(/\r?\n/), s = +L[1], cell = [2, 3, 4].map(i => L[i].trim().split(/\s+/).map(v => +v * s)), els = L[5].trim().split(/\s+/), cnt = L[6].trim().split(/\s+/).map(Number), n = cnt.reduce((a, b) => a + b, 0);
  if (!Number.isFinite(s) || cell.flat().some(v => !Number.isFinite(v)) || !(n > 0) || els.some(isNum)) throw new Error('XDATCAR: expected title, scale, three lattice vectors, element symbols and counts');
  const sym = els.flatMap((e, i) => new Array(cnt[i]).fill(e)), frames = []; let total = 0;
  for (let i = 7; i < L.length; i++) if (/^\s*Direct configuration/i.test(L[i])) { total++; if (frames.length >= 40) continue; const atoms = []; for (let k = 0; k < n; k++) { const u = (L[i + 1 + k] || '').trim().split(/\s+/).map(Number); if (u.length < 3 || u.some(v => !Number.isFinite(v))) throw new Error(`XDATCAR: configuration ${total} is truncated`); atoms.push({ el: sym[k], x: u[0] * cell[0][0] + u[1] * cell[1][0] + u[2] * cell[2][0], y: u[0] * cell[0][1] + u[1] * cell[1][1] + u[2] * cell[2][1], z: u[0] * cell[0][2] + u[1] * cell[1][2] + u[2] * cell[2][2], f: u.slice(0, 3) }); } frames.push({ t: frames.length, atoms, cell }); }
  if (!frames.length) throw new Error('XDATCAR: no "Direct configuration=" blocks found'); return { atoms: frames[0].atoms, cell, title: L[0].trim(), frames, nFrames: total, warnings: total > frames.length ? [`trajectory holds ${total} configurations; the first ${frames.length} are kept as animation frames`] : [] };
}
function strings(b, minLen = 6, limit = 10) { const out = []; let cur = ''; for (let i = 0, e = min(b.length, 65536); i < e && out.length < limit; i++) { const c = b[i]; if (c >= 32 && c < 127) cur += String.fromCharCode(c); else { if (cur.length >= minLen && /[A-Za-z]{3}/.test(cur)) out.push(cur.slice(0, 80)); cur = ''; } } return out; }
const count = (t, re) => (t.match(re) || []).length;
// Cheap, format-specific header facts for files that are recognised but not reconstructed.
function sniff(id, c) {
  const b = c.bytes, m = { bytes: b.length }, head = latin(b, 0, 4096);
  try {
    if (id === 'dwg') { const v = head.slice(0, 6); m.signature = v; m.release = ({ AC1009: 'R11/R12', AC1012: 'R13', AC1014: 'R14', AC1015: '2000', AC1018: '2004', AC1021: '2007', AC1024: '2010', AC1027: '2013', AC1032: '2018' })[v] || 'unknown release'; }
    else if (id === 'jt') m.version = (/^Version ([\d.]+) JT/.exec(head) || [])[1];
    else if (id === 'parasolid') { m.encoding = head.startsWith('PS') ? 'binary' : 'text'; for (const k of ['SCH', 'APPL', 'FORMAT', 'KEY', 'DATE', 'USFLD_SIZE']) { const r = new RegExp(k + '=([^;\\r\\n]*)').exec(head); if (r) m[k.toLowerCase()] = r[1].trim(); } }
    else if (id === 'sat') { const L = head.split(/\r?\n/), a = L[0].trim().split(/\s+/); m.version = +a[0] / 100; m.records = +a[1]; m.bodies = +a[2]; const p = /^\d+ (.*?) \d+ (ACIS [^\d]*[\d.]+[^\d]*) \d+ (.*)$/.exec(L[1] || ''); if (p) { m.product = p[1]; m.acis = p[2].trim(); m.date = p[3]; } const u = (L[2] || '').trim().split(/\s+/).map(Number); if (u.length >= 2) { m.unitScale_mm = u[0]; m.resabs = u[1]; } const t = txt(c); m.entityCounts = { body: count(t, /^-?\d*\s*body\b/gm), lump: count(t, /^-?\d*\s*lump\b/gm), shell: count(t, /^-?\d*\s*shell\b/gm), face: count(t, /^-?\d*\s*face\b/gm), edge: count(t, /^-?\d*\s*edge\b/gm), 'spline surface': count(t, /spline-surface/g), 'plane surface': count(t, /plane-surface/g), 'cone surface': count(t, /cone-surface/g) }; }
    else if (id === 'ifc') { const s = F.parseSTEP(b); Object.assign(m, s.meta); m.buildingElements = Object.fromEntries(Object.entries(s.meta.entityCounts || {}).filter(([k]) => /^IFC(WALL|SLAB|BEAM|COLUMN|PIPE|FLOW|TANK|DOOR|WINDOW|SPACE|BUILDING|EQUIPMENT|FURNISHING|DUCT|VALVE|PUMP)/.test(k))); }
    else if (id === 'brep') { m.version = (/CASCADE Topology V(\d+)/.exec(head) || [])[1]; const t = txt(c); for (const k of ['Locations', 'Curve2ds', 'Curves', 'Polygon3D', 'PolygonOnTriangulations', 'Surfaces', 'Triangulations', 'TShapes']) { const r = new RegExp('^' + k + ' (\\d+)', 'm').exec(t); if (r) m[k[0].toLowerCase() + k.slice(1)] = +r[1]; } }
    else if (id === 'las' || id === 'laz') { const v = new DataView(b.buffer, b.byteOffset, b.byteLength); if (b.length >= 227) { m.version = `${b[24]}.${b[25]}`; m.systemId = latin(b, 26, 32).replace(/\0.*$/, ''); m.software = latin(b, 58, 32).replace(/\0.*$/, ''); m.pointFormat = b[104] & 0x3F; m.compressed = !!(b[104] & 0x80) || id === 'laz'; m.pointRecordBytes = v.getUint16(105, true); m.points = v.getUint32(107, true); m.scale = [131, 139, 147].map(o => v.getFloat64(o, true)); m.offset = [155, 163, 171].map(o => v.getFloat64(o, true)); m.bboxMax = [179, 195, 211].map(o => v.getFloat64(o, true)); m.bboxMin = [187, 203, 219].map(o => v.getFloat64(o, true)); if (b[24] === 1 && b[25] >= 4 && b.length >= 255 && !m.points) m.points = Number(v.getBigUint64(247, true)); } }
    else if (id === 'e57') { const v = new DataView(b.buffer, b.byteOffset, b.byteLength); m.version = `${v.getUint32(8, true)}.${v.getUint32(12, true)}`; }
    else if (id === 'xtc') { const v = new DataView(b.buffer, b.byteOffset, b.byteLength), mg = v.getInt32(0, false); m.kind = mg === 1995 ? 'XTC (compressed coordinates)' : 'TRR (full precision)'; if (mg === 1995) { m.atoms = v.getInt32(4, false); m.firstStep = v.getInt32(8, false); m.firstTime_ps = v.getFloat32(12, false); } }
    else if (id === 'vdb') { const v = new DataView(b.buffer, b.byteOffset, b.byteLength); m.fileVersion = v.getUint32(8, true); m.embeddedStrings = strings(b); }
    else if (id === 'fbx') { m.encoding = head.startsWith('Kaydara') ? 'binary' : 'ASCII'; if (m.encoding === 'binary' && b.length > 27) m.version = new DataView(b.buffer, b.byteOffset, b.byteLength).getUint32(23, true); }
    else if (id === 'usd') { m.variant = head.startsWith('PXR-USDC') ? 'crate (binary)' : head.startsWith('#usda') ? 'text' : 'package'; if (m.variant === 'text') { const t = txt(c); m.prims = { Mesh: count(t, /def Mesh\b/g), Xform: count(t, /def Xform\b/g), Scope: count(t, /def Scope\b/g), Material: count(t, /def Material\b/g) }; m.upAxis = (/upAxis\s*=\s*"(\w)"/.exec(t) || [])[1]; m.metersPerUnit = +((/metersPerUnit\s*=\s*([\d.eE+-]+)/.exec(t) || [])[1]) || null; } }
    else if (id === 'dae') { const t = txt(c); m.geometries = count(t, /<geometry\b/g); m.triangleBlocks = count(t, /<(triangles|polylist|polygons)\b/g); m.unit = (/<unit\b[^>]*meter="([^"]+)"/.exec(t) || [])[1]; m.upAxis = (/<up_axis>(\w+)</.exec(t) || [])[1]; }
    else if (id === 'fluent') { if (!head.startsWith('\x89HDF')) { const t = txt(c); m.dimension = +((/\(\s*2\s+(\d)\s*\)/.exec(t) || [])[1]) || null; const n = /\(\s*10\s+\(\s*0\s+\w+\s+(\w+)/.exec(t), f = /\(\s*13\s+\(\s*0\s+\w+\s+(\w+)/.exec(t), cl = /\(\s*12\s+\(\s*0\s+\w+\s+(\w+)/.exec(t); m.nodes = n ? parseInt(n[1], 16) : null; m.faces = f ? parseInt(f[1], 16) : null; m.cells = cl ? parseInt(cl[1], 16) : null; m.zones = [...t.matchAll(/\(\s*(?:39|45)\s+\(\s*(\d+)\s+([\w-]+)\s+([^\s)]+)/g)].slice(0, 60).map(z => ({ id: +z[1], type: z[2], name: z[3] })); } }
    else if (id === 'cdb') { const t = txt(c), nb = /NBLOCK,\s*(\d+),\s*\w*,\s*(\d+),\s*(\d+)/.exec(t), eb = /EBLOCK,\s*(\d+),\s*\w*,\s*(\d+),\s*(\d+)/.exec(t); m.nodes = nb ? +nb[3] : null; m.elements = eb ? +eb[3] : null; m.elementTypes = [...t.matchAll(/^ET,\s*(\d+),\s*(\w+)/gm)].slice(0, 30).map(e => `${e[1]}: ${e[2]}`); m.components = count(t, /^CMBLOCK,/gm); }
    else if (id === 'lsdyna') { const t = txt(c), kw = {}; for (const k of t.matchAll(/^\*([A-Z_0-9]+)/gm)) kw[k[1]] = (kw[k[1]] || 0) + 1; m.keywords = Object.fromEntries(Object.entries(kw).sort((a, b2) => b2[1] - a[1]).slice(0, 30)); m.title = (/^\*TITLE\s*\r?\n([^\r\n]*)/m.exec(t) || [])[1]; }
    else if (id === 'ensight') { const t = txt(c); m.model = (/^model:\s*(?:\d+\s+)*(\S+)/m.exec(t) || [])[1]; m.variables = [...t.matchAll(/^(scalar|vector|tensor \w+) per (node|element):\s*(?:\d+\s+)*(\S+)\s+(\S+)/gm)].slice(0, 30).map(v => `${v[3]} (${v[1]} per ${v[2]}) ← ${v[4]}`); m.timeSteps = +((/number of steps:\s*(\d+)/.exec(t) || [])[1]) || null; }
    else if (id === 'openfoam') { const t = txt(c); m.class = (/class\s+(\w+);/.exec(t) || [])[1]; m.object = (/object\s+(\w+);/.exec(t) || [])[1]; m.entries = +((/^\s*(\d+)\s*\r?\n\s*\(/m.exec(t) || [])[1]) || null; if (m.object === 'boundary') m.patches = [...t.matchAll(/^\s{4}(\w+)\s*\r?\n\s*\{[^}]*type\s+(\w+);/gm)].map(p => `${p[1]} (${p[2]})`); }
    else if (id === 'lef' || id === 'def') { const t = txt(c); if (id === 'lef') { m.macros = [...t.matchAll(/^\s*MACRO\s+(\S+)/gm)].slice(0, 100).map(x => x[1]); m.layers = [...new Set([...t.matchAll(/^\s*LAYER\s+(\S+)\s*$/gm)].map(x => x[1]))].slice(0, 60); m.databaseMicrons = +((/DATABASE MICRONS\s+(\d+)/.exec(t) || [])[1]) || null; } else { m.design = (/^\s*DESIGN\s+(\S+)/m.exec(t) || [])[1]; m.unitsPerMicron = +((/UNITS DISTANCE MICRONS\s+(\d+)/.exec(t) || [])[1]) || null; m.dieArea = (/DIEAREA\s+([^;]*);/.exec(t) || [])[1]; m.components = +((/^\s*COMPONENTS\s+(\d+)/m.exec(t) || [])[1]) || 0; m.nets = +((/^\s*NETS\s+(\d+)/m.exec(t) || [])[1]) || 0; m.pins = +((/^\s*PINS\s+(\d+)/m.exec(t) || [])[1]) || 0; } }
    else if (id === 'caltechcif') { const t = txt(c); m.symbolDefinitions = count(t, /^\s*DS\s+\d+/gm); m.boxes = count(t, /(^|;)\s*B\s+-?\d/gm); m.polygons = count(t, /(^|;)\s*P\s+-?\d/gm); m.wires = count(t, /(^|;)\s*W\s+-?\d/gm); m.calls = count(t, /(^|;)\s*C\s+\d/gm); m.layers = [...new Set([...t.matchAll(/(?:^|;)\s*L\s+(\w+)/gm)].map(x => x[1]))]; m.units = 'centimicrons (0.01 µm) scaled by the DS a/b factors'; }
    else if (id === 'python') { const t = txt(c); m.lines = t.split('\n').length; m.imports = [...new Set([...t.matchAll(/^\s*(?:import|from)\s+([\w.]+)/gm)].map(x => x[1].split('.')[0]))].slice(0, 40); m.functions = [...t.matchAll(/^\s*def\s+(\w+)\s*\(([^)]*)\)/gm)].slice(0, 60).map(x => `${x[1]}(${x[2].replace(/\s+/g, ' ').slice(0, 60)})`); m.geometryLibraries = m.imports.filter(i => /^(cadquery|build123d|gdspy|gdstk|gdsfactory|ase|pymatgen|trimesh|pygmsh|gmsh|numpy|scipy|FreeCAD|Part|bpy|nazca|phidl|meep|shapely|solid|openpyscad|pyvista|skimage|porespy|meshio)$/.test(i)); m.executed = false; }
    else if (id === 'grasshopper') { if (/<Archive/.test(head)) { const t = txt(c); m.encoding = 'ghx (XML)'; m.objects = +((/<item name="ObjectCount"[^>]*>(\d+)</.exec(t) || [])[1]) || count(t, /<chunk name="Object"/g); m.components = Object.fromEntries(Object.entries([...t.matchAll(/<item name="Name" type_name="gh_string"[^>]*>([^<]+)</g)].reduce((o, x) => { o[x[1]] = (o[x[1]] || 0) + 1; return o; }, {})).sort((a, b2) => b2[1] - a[1]).slice(0, 30)); } else { m.encoding = 'gh (binary archive)'; m.embeddedStrings = strings(b); } }
    else if (id === 'dynamo') { const t = txt(c); try { const j = JSON.parse(t); m.name = j.Name; m.nodes = (j.Nodes || []).length; m.connectors = (j.Connectors || []).length; m.nodeTypes = Object.entries((j.Nodes || []).reduce((o, n) => { const k = n.ConcreteType ? n.ConcreteType.split(',')[0].split('.').pop() : n.NodeType || 'node'; o[k] = (o[k] || 0) + 1; return o; }, {})).slice(0, 30).map(([k, v]) => `${k} × ${v}`); } catch (e) { m.nodes = count(t, /<Dynamo\.[\w.]+|<[\w.]*Node\b/g); } }
    else if (id === 'zarr' || id === 'omezarr') { try { const j = JSON.parse(txt(c)); Object.assign(m, { zarrFormat: j.zarr_format, shape: j.shape, chunks: j.chunks, dtype: j.dtype, compressor: j.compressor && j.compressor.id, order: j.order }); } catch (e) { m.note = 'archive or directory store'; } }
    else if (id === 'mmcif') { const t = txt(c); m.entry = (/^_entry\.id\s+(\S+)/m.exec(t) || /^data_(\S+)/m.exec(t) || [])[1]; m.atomRecords = count(t, /^(ATOM|HETATM)\s/gm); m.cell = ['length_a', 'length_b', 'length_c', 'angle_alpha', 'angle_beta', 'angle_gamma'].map(k => +((new RegExp('_cell\\.' + k + '\\s+([\\d.]+)').exec(t) || [])[1])); }
    else if (id === 'qeout') { const t = txt(c); m.program = (/Program (PWSCF[^\r\n]*)/.exec(t) || [])[1]; m.atoms = +((/number of atoms\/cell\s*=\s*(\d+)/.exec(t) || [])[1]) || null; const en = [...t.matchAll(/^!\s+total energy\s+=\s+([-\d.]+) Ry/gm)]; m.finalEnergy_Ry = en.length ? +en[en.length - 1][1] : null; m.relaxationSteps = en.length; }
    else if (id === 'abinit') { const t = txt(c); for (const k of ['natom', 'ntypat', 'acell', 'ecut', 'znucl']) { const r = new RegExp('^\\s*' + k + '\\s+([^\\r\\n#]*)', 'm').exec(t); if (r) m[k] = r[1].trim(); } }
    else if (id === 'shelx') { const t = txt(c); m.title = (/^TITL\s+(.*)$/m.exec(t) || [])[1]; m.cell = (/^CELL\s+(.*)$/m.exec(t) || [])[1]; m.elements = (/^SFAC\s+(.*)$/m.exec(t) || [])[1]; m.latt = (/^LATT\s+(.*)$/m.exec(t) || [])[1]; m.symmetryOperators = count(t, /^SYMM\s/gm); }
    else if (id === 'psf' || id === 'gmxtop' || id === 'charmmcrd' || id === 'zmat' || id === 'wien2k') { const t = txt(c); m.headerLines = t.split(/\r?\n/).filter(l => l.trim()).slice(0, 4).map(l => l.slice(0, 120)); if (id === 'psf') { m.atoms = +((/^\s*(\d+)\s+!NATOM/m.exec(t) || [])[1]) || null; m.bonds = +((/^\s*(\d+)\s+!NBOND/m.exec(t) || [])[1]) || null; } if (id === 'gmxtop') m.moleculeTypes = [...t.matchAll(/\[\s*moleculetype\s*\]\s*(?:;[^\n]*\n|\s)*(\S+)/g)].map(x => x[1]).slice(0, 40); }
    else if (id === 'plot3d') { const a = head.trim().split(/\s+/); if (a.slice(0, 4).every(s => /^\d+$/.test(s))) { m.blocks = +a[0]; m.firstBlockDims = a.slice(1, 4).map(Number); } else if (a.slice(0, 3).every(s => /^\d+$/.test(s))) m.dims = a.slice(0, 3).map(Number); else m.encoding = 'binary (Fortran unformatted or C stream)'; }
    else if (id === 'slc') { const e = head.indexOf('\x1a'); m.header = head.slice(0, e > 0 ? e : 300).replace(/[\r\n]+/g, ' ').trim().slice(0, 300); m.units = (/-UNIT\s+(\w+)/.exec(head) || [])[1]; m.extents = (/-EXTENTS\s+([^-]*)/.exec(head) || [])[1]; }
    else if (id === 'vtm') { const t = txt(c); m.blocks = [...t.matchAll(/<DataSet\b([^>]*)\/?>/g)].map(x => (/file="([^"]*)"/.exec(x[1]) || [])[1]).filter(Boolean).slice(0, 100); }
    else if (id === 'mtex') { const t = F.parseTable(txt(c)); m.columns = t.header; m.rows = t.rows.length; }
    else if (id === 'xml') { const t = txt(c); m.root = (/<([A-Za-z][\w:.-]*)[\s>]/.exec(t.replace(/<\?[\s\S]*?\?>|<!--[\s\S]*?-->|<![^>]*>/g, '')) || [])[1]; const tags = {}; for (const x of t.matchAll(/<([A-Za-z][\w:.-]*)[\s>\/]/g)) tags[x[1]] = (tags[x[1]] || 0) + 1; m.tagCounts = Object.fromEntries(Object.entries(tags).sort((a, b2) => b2[1] - a[1]).slice(0, 25)); }
    else if (id === 'fcstd' || id === 'fusion' || id === 'odbpp' || (head.startsWith('PK') && ['solidworks', 'usd', 'unknown'].includes(id))) { const en = F.unzip(b, /^Document\.xml$/); m.archiveEntries = en.slice(0, 40).map(e => `${e.name} (${e.size} bytes)`); const d = en.find(e => e.name === 'Document.xml' && e.data); if (d) { const t = F.toText(d.data); m.programVersion = (/ProgramVersion="([^"]+)"/.exec(t) || [])[1]; m.objects = [...t.matchAll(/<Object type="([^"]+)" name="([^"]+)"/g)].slice(0, 80).map(o => `${o[2]} : ${o[1]}`); } }
    else if (id === 'dicom') { m.transferSyntaxHint = strings(b.subarray(128, 1024), 8, 6); }
    else if (['cgns', 'med', 'h5ebsd', 'exodus', 'abc', 'omezarr'].includes(id) && b[0] === 0x89) { const h = F.h5Open(b); m.groups = h.groups.slice(0, 60); m.datasets = h.datasets.slice(0, 80).map(d => `${d.path} ${d.dtype} [${d.shape.join(' × ')}]`); m.nDatasets = h.datasets.length; if (h.notes.length) m.notes = h.notes; }
    else if (!Object.keys(m).some(k => k !== 'bytes')) { let nul = false; for (let i = 0; i < min(b.length, 2048); i++) if (b[i] === 0) nul = true; if (nul) { m.magicHex = Array.from(b.subarray(0, 16)).map(v => v.toString(16).padStart(2, '0')).join(' '); m.embeddedStrings = strings(b); } else m.headerLines = head.split(/\r?\n/).filter(l => l.trim()).slice(0, 5).map(l => l.slice(0, 120)); }
  } catch (e) { m.sniffError = e.message; }
  return m;
}

// ---------------------------------------------------------------- handlers: id → geometry (or null) with notes
export const GENERATORS = {
  tpms: p => tpmsVoxels(p), lattice: p => buildLattice(p), porous: p => porousMedium(p), porousmedium: p => porousMedium(p), nanorods: p => nanorodArray(p), nanorodarray: p => nanorodArray(p), reactor: p => { const r = reactor(p); return { ...r.mesh, regions: r.regions, volume_m3: r.volume_m3, wettedArea_m2: r.wettedArea_m2 }; },
  resonator: p => layoutArray(resonator(p), { nx: 1, ny: 1, ...(p.array || {}) }), metasurface: p => layoutArray(resonator(p), p.array || {}),
};
const struct = (obj, what) => { const g = F.fromStructured(obj, GENERATORS); g.notes = [...(g.notes || [])]; g.meta = { ...(g.meta || {}), source: what }; return g; };
const atomic = fn => c => atomsOut(fn(txt(c), c.opts));
const H = {
  stl: c => F.parseSTL(c.bytes), obj: c => F.parseOBJ(c.bytes), ply: c => F.parsePLY(c.bytes), off: c => F.parseOFF(c.bytes), gts: c => F.parseGTS(c.bytes), '3mf': c => F.parse3MF(c.bytes), amf: c => F.parseAMF(c.bytes),
  gltf: c => F.parseGLTF(c.bytes), vrml: c => F.parseVRML(c.bytes), x3d: c => F.parseVRML(c.bytes), step: c => F.parseSTEP(c.bytes), iges: c => F.parseIGES(c.bytes),
  vtk: c => F.parseVTK(c.bytes), vtu: c => F.parseVTU(c.bytes), vtp: c => F.parseVTP(c.bytes), vti: c => F.parseVTI(c.bytes), msh: c => F.parseMSH(c.bytes), inp: c => F.parseINP(c.bytes), bdf: c => F.parseBDF(c.bytes), unv: c => F.parseUNV(c.bytes), su2: c => F.parseSU2(c.bytes),
  xdmf: c => F.parseXDMF(c.bytes, { files: c.opts.files || {} }), pvd: c => F.parsePVD(c.bytes, { files: c.opts.files || {} }),
  gds: c => F.parseGDS(c.bytes), oasis: c => F.parseOASIS(c.bytes), dxf: c => F.parseDXF(c.bytes), svg: c => F.parseSVG(c.bytes), gerber: c => F.parseGerber(c.bytes), excellon: c => F.parseExcellon(c.bytes),
  cli: c => F.parseCLI(c.bytes), gcode: c => F.parseGcode(c.bytes), scad: c => F.parseSCAD(c.bytes, { resolution: c.opts.resolution || 64 }), geo: c => F.parseGEO(c.bytes),
  tiff: c => F.parseTIFF(c.bytes, c.opts), ometiff: c => F.parseTIFF(c.bytes, c.opts), hdf5: c => F.parseHDF5(c.bytes, c.opts), nrrd: c => F.parseNRRD(c.bytes, c.opts), metaimage: c => F.parseMHD(c.bytes, c.opts), nifti: c => F.parseNIfTI(c.bytes, c.opts), mrc: c => F.parseMRC(c.bytes, c.opts), raw: c => F.parseRAW(c.bytes, c.opts),
  em: c => { const b = c.bytes; if (b.length < 512) throw new Error('EM: file shorter than the 512-byte header'); const v = new DataView(b.buffer, b.byteOffset, b.byteLength), le = b[0] === 1 || b[0] === 6 /* machine code: 1 VAX, 6 PC are little-endian; 0, 2, 3, 5 big-endian */, dt = ({ 1: 'int8', 2: 'int16', 4: 'int32', 5: 'float32', 9: 'float64' })[b[3]]; if (!dt) throw new Error(`EM: data-type code ${b[3]} is not supported`); const g = F.parseRAW(b, { nx: v.getInt32(4, le), ny: v.getInt32(8, le), nz: v.getInt32(12, le), dtype: dt, littleEndian: le, headerBytes: 512, threshold: c.opts.threshold }); g.notes = g.notes.filter(n => !/raw volumes carry no metadata/.test(n)).concat('EM header gives dimensions and sample type only: voxel size is not calibrated'); return g; },
  ctf: c => F.parseCTF(c.bytes), ang: c => F.parseANG(c.bytes), graphml: c => F.parseGraphML(c.bytes), gexf: c => F.parseGEXF(c.bytes), gml: c => F.parseGML(c.bytes),
  xyzpoints: c => F.parsePointCloud(c.bytes), pts: c => F.parsePointCloud(c.bytes), ptx: c => F.parsePointCloud(c.bytes), pcd: c => F.parsePointCloud(c.bytes),
  json: c => { let j; try { j = JSON.parse(txt(c)); } catch (e) { throw new Error('JSON: ' + e.message); } return struct(j, 'JSON'); }, yaml: c => struct(F.parseYAML(txt(c)), 'YAML'),
  csv: c => { const t = txt(c), tb = F.parseTable(t), h = tb.header || []; if (h.some(k => /thick/.test(k))) { const mi = h.findIndex(k => /mat|name|layer/.test(k)), ti = h.findIndex(k => /thick/.test(k)), u = (/\((\w+)\)|_(nm|um|mm)\b/.exec(h[ti]) || []).slice(1).find(Boolean); return struct({ units: u || 'nm', layers: tb.rows.map(r => ({ material: mi >= 0 ? r[mi] : undefined, thickness: +r[ti] })) }, 'CSV layer table'); } if (h.some(k => ['source', 'from', 'target', 'x1', 'node1', 'n1'].includes(k)) || t.split(/\r?\n/).some(l => F.EDGE_HEADER.test(l))) return F.parseNodeEdgeCSV(t); return F.parsePointCloud(c.bytes); },
  xml: c => { const t = txt(c), ls = [...t.matchAll(/<layer\b([^>]*)\/?>/gi)].map(m => ({ material: (/\b(?:material|name)\s*=\s*"([^"]*)"/i.exec(m[1]) || [])[1], thickness: +((/\bthickness(?:_nm)?\s*=\s*"([^"]*)"/i.exec(m[1]) || [])[1]) })).filter(l => Number.isFinite(l.thickness)); if (ls.length) return struct({ layers: ls, units: (/\bunits?\s*=\s*"(\w+)"/i.exec(t) || [])[1] }, 'XML layer table'); return null; },
  xyz: atomic(A.parseXYZ), extxyz: atomic(A.parseXYZ), pdb: atomic(A.parsePDB), pqr: atomic(A.parsePDB), mol: atomic(A.parseMOL), mol2: atomic(A.parseMOL2), cml: atomic(A.parseCML), gaussian: atomic(A.parseGaussianInput), gro: atomic(A.parseGRO),
  sdf: c => { const t = txt(c), n = (t.match(/^\$\$\$\$/gm) || []).length, g = atomsOut(A.parseSDF(t)); if (n > 1) g.lost.push(`${n - 1} further molecule records in the SD file were not imported`); return g; },
  lammpsdata: atomic(A.parseLAMMPSData), lammpsdump: c => atomsOut(A.parseLAMMPSDump(txt(c), { frames: 'all', maxFrames: 40, ...c.opts })), poscar: atomic(A.parsePOSCAR), xsf: atomic(A.parseXSF), cif: atomic(A.parseCIF), qein: atomic(A.parseQEInput), castep: atomic(A.parseCASTEPCell),
  fhiaims: c => atomsOut(parseAims(txt(c))), cp2k: c => atomsOut(parseCP2K(txt(c)), [], ['basis sets, functionals and run settings are not geometry and were ignored']), xdatcar: c => atomsOut(parseXDATCAR(txt(c))),
  dream3d: c => { const h = F.h5Open(c.bytes), d = h.datasets.find(q => q.readable && /\/Phases$/.test(q.path) && q.shape.filter(n => n > 1).length >= 2) || h.datasets.find(q => q.readable && /uint8|int8/.test(q.dtype) && q.shape.filter(n => n > 1).length === 3); if (!d) return null; const g = F.parseHDF5(c.bytes, { ...c.opts, dataset: d.path }); g.lost.push('FeatureIds, Euler angles / quaternions, grain statistics and pipeline metadata were not imported'); return g; },
};

// ---------------------------------------------------------------- validation, periodicity, suitability
const chk = (name, pass, detail, info) => ({ name, pass: !!pass, detail, ...(info ? { info: true } : {}) }); // info: a property of the geometry, not a defect
// Translational periodicity of a point set: points on the low face of an axis map onto points on the high face.
export function detectPeriodicity(points, maxPoints = 200000) {
  if (!points || points.length < 4 || points.length > maxPoints) return null;
  const [lo, hi] = bboxOf(points), size = [0, 1, 2].map(a => hi[a] - lo[a]), diag = hypot(...size) || 1, tol = 1e-6 * diag, key = (a, p) => { const b = (a + 1) % 3, c = (a + 2) % 3; return round((p[b] || 0) / (10 * tol)) + ',' + round((p[c] || 0) / (10 * tol)); }, out = { axes: {}, period: size, repeats: {} };
  for (let a = 0; a < 3; a++) {
    if (size[a] < tol) { out.axes['xyz'[a]] = false; continue; }
    const L = new Set(), Hs = new Set(); for (const p of points) { const v = p[a] || 0; if (v - lo[a] < tol) L.add(key(a, p)); else if (hi[a] - v < tol) Hs.add(key(a, p)); }
    let match = 0; for (const k of L) if (Hs.has(k)) match++;
    out.axes['xyz'[a]] = L.size > 0 && L.size === Hs.size && match === L.size;
    // repeated sub-cells: does a shift by size/n map the point set onto itself (points near the far end excepted)?
    if (points.length <= 60000) { const all = new Set(points.map(p => round((p[0] || 0) / (10 * tol)) + ',' + round((p[1] || 0) / (10 * tol)) + ',' + round((p[2] || 0) / (10 * tol)))); for (let n = 8; n >= 2; n--) { const s = size[a] / n; let ok = 0, tot = 0; for (const p of points) { if ((p[a] || 0) + s > hi[a] + tol) continue; tot++; const q = [p[0] || 0, p[1] || 0, p[2] || 0]; q[a] += s; if (all.has(round(q[0] / (10 * tol)) + ',' + round(q[1] / (10 * tol)) + ',' + round(q[2] / (10 * tol)))) ok++; } if (tot > 0 && ok === tot) { out.repeats['xyz'[a]] = n; break; } } }
  }
  out.any = Object.values(out.axes).some(Boolean); return out;
}
function validate(g, rep, extra, opts) {
  const V = rep.validation, facts = {};
  if (g.kind === 'mesh') {
    if (!g.vertices.length || !g.faces.length) { V.push(chk('Mesh content', false, 'no vertices or faces')); return facts; }
    const bad = g.faces.some(f => f.length !== 3 || f.some(i => !(i >= 0 && i < g.vertices.length)) ) || g.vertices.some(p => !Number.isFinite(p[0] + p[1] + (p[2] || 0)));
    V.push(chk('Index and coordinate integrity', !bad, bad ? 'faces reference missing vertices or coordinates are not finite' : `${g.vertices.length} vertices, ${g.faces.length} triangles, all indices valid`)); if (bad) return facts;
    let st = meshStats(g); const before = st, shells = !!g.shells;
    if (opts.repair !== false && !shells && g.faces.length <= 1.5e6 && (st.duplicateVertices || st.degenerateFaces || !st.consistentNormals || st.outwardNormals === false || st.unreferencedVertices)) { const r = repairMesh(g); const changed = r.report.filter(s => !/^no defects|remain/.test(s)); if (changed.length) { g.vertices = r.mesh.vertices; g.faces = r.mesh.faces; rep.repairs.push(...changed); st = r.stats; } }
    V.push(chk('Duplicate vertices', st.duplicateVertices === 0, st.duplicateVertices ? `${st.duplicateVertices} coincident vertices${shells ? ' (expected: overlapping shells)' : ''}` : before.duplicateVertices ? `${before.duplicateVertices} found and welded` : 'none'));
    V.push(chk('Degenerate faces', st.degenerateFaces === 0, st.degenerateFaces ? `${st.degenerateFaces} zero-area faces` : before.degenerateFaces ? `${before.degenerateFaces} found and removed` : 'none'));
    V.push(chk('Open edges (watertight)', st.boundaryEdges === 0, st.boundaryEdges ? `${st.boundaryEdges} boundary edges: the surface is open` : 'every edge is shared by two faces'));
    V.push(chk('Non-manifold edges', st.nonManifoldEdges === 0, st.nonManifoldEdges ? `${st.nonManifoldEdges} edges shared by more than two faces` : 'none'));
    V.push(chk('Consistent normals', st.consistentNormals, st.consistentNormals ? (before.consistentNormals ? 'winding consistent' : 'winding made consistent by repair') : `${st.inconsistentEdges} edges traversed twice in the same direction`));
    if (st.watertight) V.push(chk('Outward orientation', st.volume > 0, `signed volume ${st.volume.toPrecision(5)} ${g.units || 'units'}³`));
    rep.assumptions.push('self-intersections are not tested');
    Object.assign(facts, { watertight: st.watertight, manifold: st.nonManifoldEdges === 0, volume: st.volume, solid: !!g.solidElements, components: st.components }); extra.stats = st;
    const per = detectPeriodicity(g.vertices); if (per) { extra.periodicity = per; facts.periodic = per.any; V.push(chk('Periodic boundary match', per.any, per.any ? `matching vertex sets on opposite faces along ${Object.entries(per.axes).filter(e => e[1]).map(e => e[0]).join(', ')}` : 'opposite bounding-box faces do not carry matching vertices (not a periodic unit cell as meshed)', true)); }
  } else if (g.kind === 'atoms') {
    const n = g.atoms.length; V.push(chk('Atom records', n > 0 && g.atoms.every(a => a.el && Number.isFinite(a.x + a.y + a.z)), `${n} atoms with element symbol and finite coordinates`)); if (!n) return facts;
    if (g.cell) { const c = g.cell, vol = abs(c[0][0] * (c[1][1] * c[2][2] - c[1][2] * c[2][1]) - c[0][1] * (c[1][0] * c[2][2] - c[1][2] * c[2][0]) + c[0][2] * (c[1][0] * c[2][1] - c[1][1] * c[2][0])); V.push(chk('Lattice vectors', vol > 1e-6, `cell volume ${vol.toPrecision(6)} Å³`)); facts.periodic = vol > 1e-6; V.push(chk('Periodicity', true, 'periodic cell present: lattice vectors preserved')); } else V.push(chk('Periodicity', true, 'no cell: finite (non-periodic) structure'));
    if (n <= 4000) { let dmin = Infinity; const at = g.atoms; for (let i = 0; i < n; i++) for (let j = i + 1; j < n; j++) { const d = (at[i].x - at[j].x) ** 2 + (at[i].y - at[j].y) ** 2 + (at[i].z - at[j].z) ** 2; if (d < dmin) dmin = d; } dmin = Math.sqrt(dmin); if (n > 1) V.push(chk('Overlapping atoms', dmin > 0.5, `shortest interatomic distance ${dmin.toFixed(3)} Å (minimum-image contacts across the cell are not included)`)); } else rep.assumptions.push('interatomic-distance check skipped above 4000 atoms');
    try { if (n <= 6000) { const an = A.analyse(g); extra.stats = { formula: an.formula, composition: an.composition, density_kg_m3: an.density_kg_m3, volume_A3: an.volume_A3, nearestNeighbour_A: an.nearestNeighbour_A }; rep.imported.push(`composition ${an.formula}${an.density_kg_m3 ? `, density ${Math.round(an.density_kg_m3)} kg/m³` : ''}`); } } catch (e) { /* analysis is optional */ }
    V.push(chk('Fractional coordinates', true, g.atoms[0].f ? 'fractional coordinates kept alongside Cartesian positions' : 'Cartesian positions only (fractional coordinates can be recomputed from the cell)'));
    facts.atoms = n;
  } else if (g.kind === 'voxels') {
    const N = g.nx * g.ny * g.nz, ok = g.data && g.data.length === N && N > 0; V.push(chk('Grid consistency', ok, ok ? `${g.nx} × ${g.ny} × ${g.nz} = ${N} voxels` : `data length ${g.data && g.data.length} ≠ nx·ny·nz = ${N}`)); if (!ok) return facts;
    if (N <= 4e6) { const s = voxelStats(g); extra.stats = s; const ph = Object.keys(s.fractions); V.push(chk('Phase content', ph.length > 1, ph.length > 1 ? `${ph.length} phase ids; fractions ${ph.slice(0, 6).map(k => `${k}: ${(s.fractions[k] * 100).toFixed(1)} %`).join(', ')}` : `single phase id ${ph[0]}: the volume is uniform (check the threshold)`)); if (s.percolates) { V.push(chk('Pore-phase percolation (phase 0)', s.percolates.x || s.percolates.y || s.percolates.z, `x: ${s.percolates.x}, y: ${s.percolates.y}, z: ${s.percolates.z}`, true)); V.push(chk('Solid-phase connectivity', s.solidPercolates.x || s.solidPercolates.y || s.solidPercolates.z, `x: ${s.solidPercolates.x}, y: ${s.solidPercolates.y}, z: ${s.solidPercolates.z}`, true)); Object.assign(facts, { porePercolates: s.percolates.x || s.percolates.y || s.percolates.z, solidPercolates: s.solidPercolates.x || s.solidPercolates.y || s.solidPercolates.z }); } facts.phases = ph.length; } else rep.assumptions.push('phase and percolation statistics skipped above 4 million voxels');
    V.push(chk('Voxel calibration', !!g.units && g.units !== 'px', g.units && g.units !== 'px' ? `voxel size ${g.voxel} ${g.units}` : 'no physical voxel size: lengths are in voxels until a scale is set')); facts.calibrated = !!g.units && g.units !== 'px'; facts.thin = g.nz === 1;
  } else if (g.kind === 'lattice') {
    const n = g.nodes.length, bad = g.edges.some(e => !(e[0] >= 0 && e[0] < n && e[1] >= 0 && e[1] < n)); V.push(chk('Connectivity integrity', !bad && n > 0, bad ? 'edges reference missing nodes' : `${n} nodes, ${g.edges.length} edges`)); if (bad || !n) return facts;
    if (g.pointCloud || !g.edges.length) { facts.pointCloud = true; V.push(chk('Surface information', false, 'point set without connectivity: surface reconstruction is required before solid or mesh-based analysis')); return facts; }
    const deg = new Int32Array(n), par = new Int32Array(n).map((_, i) => i), find = i => { while (par[i] !== i) { par[i] = par[par[i]]; i = par[i]; } return i; }; let zero = 0; const seen = new Set(); let dupE = 0;
    for (const e of g.edges) { deg[e[0]]++; deg[e[1]]++; const a = g.nodes[e[0]], b = g.nodes[e[1]]; if (hypot(a[0] - b[0], a[1] - b[1], (a[2] || 0) - (b[2] || 0)) < 1e-12) zero++; const k = e[0] < e[1] ? e[0] * n + e[1] : e[1] * n + e[0]; if (seen.has(k)) dupE++; else seen.add(k); const ra = find(e[0]), rb = find(e[1]); if (ra !== rb) par[ra] = rb; }
    let dangling = 0, isolated = 0; const comps = new Set(); for (let i = 0; i < n; i++) { if (deg[i] === 1) dangling++; if (deg[i] === 0) isolated++; else comps.add(find(i)); }
    const ks = new Set(); let dupN = 0; const [lo, hi] = bboxOf(g.nodes), q = 1e-7 * (hypot(hi[0] - lo[0], hi[1] - lo[1], hi[2] - lo[2]) || 1); for (const p of g.nodes) { const k = round(p[0] / q) + ',' + round(p[1] / q) + ',' + round((p[2] || 0) / q); if (ks.has(k)) dupN++; else ks.add(k); }
    V.push(chk('Duplicate nodes', dupN === 0, dupN ? `${dupN} coincident nodes (struts meeting there are not connected in the graph)` : 'none'));
    V.push(chk('Duplicate or zero-length struts', zero + dupE === 0, zero + dupE ? `${dupE} repeated and ${zero} zero-length edges` : 'none'));
    V.push(chk('Connected network', comps.size === 1, `${comps.size} connected component(s); ${dangling} dangling strut ends; ${isolated} unconnected nodes`));
    V.push(chk('Strut section defined', (g.radius || 0) > 0 || !!g.radii, (g.radius || 0) > 0 || g.radii ? `radius ${g.radii ? 'per strut' : g.radius} ${g.units || ''}` : 'no strut radius in the file: a cross-section must be assigned before mechanical analysis'));
    const Z = 2 * g.edges.length / max(1, n - isolated), dim = g.nodes.every(p => abs(p[2] || 0) < 1e-12) ? 2 : 3, M = dim === 3 ? g.edges.length - 3 * (n - isolated) + 6 : g.edges.length - 2 * (n - isolated) + 3;
    extra.stats = { nodes: n, edges: g.edges.length, components: comps.size, danglingEnds: dangling, meanConnectivity: Z, dimension: dim, maxwell: M }; Object.assign(facts, { connected: comps.size === 1, hasSection: (g.radius || 0) > 0 || !!g.radii, wireframe: !!g.wireframe, toolpath: !!g.toolpath });
    const per = detectPeriodicity(g.nodes); if (per) { extra.periodicity = per; facts.periodic = per.any; V.push(chk('Periodic unit cell', per.any, per.any ? `matching nodes on opposite faces along ${Object.entries(per.axes).filter(e => e[1]).map(e => e[0]).join(', ')}${Object.keys(per.repeats).length ? `; repeated ${Object.entries(per.repeats).map(([k, v]) => `${v}× along ${k}`).join(', ')}` : ''}` : 'opposite faces do not carry matching nodes', true)); }
  } else if (g.kind === 'layout') {
    const P = g.polygons, bad = P.filter(q => !q.pts || q.pts.length < 3 || q.pts.some(p => !Number.isFinite(p[0] + p[1]))).length, zero = P.filter(q => q.pts && q.pts.length >= 3 && abs(polyArea(q.pts)) < 1e-12).length;
    V.push(chk('Polygon integrity', P.length > 0 && bad === 0, bad ? `${bad} polygons with fewer than three finite vertices` : `${P.length} polygons`)); if (!P.length || bad) return facts;
    V.push(chk('Zero-area polygons', zero === 0, zero ? `${zero} degenerate polygons` : 'none'));
    const s = layoutStats(g); extra.stats = s; extra.layers = Object.entries(s.layers).map(([k, v]) => ({ layer: +k, polygons: v.polygons, area: v.area, name: g.layerNames ? g.layerNames[k] : undefined }));
    V.push(chk('Layer assignment', true, `${Object.keys(s.layers).length} layer(s): ${Object.entries(s.layers).slice(0, 8).map(([k, v]) => `${k} (${v.polygons})`).join(', ')}`));
    V.push(chk('Units', !!g.units && !/unitless|px/.test(g.units), g.units && !/unitless|px/.test(g.units) ? `coordinates in ${g.units}; smallest feature ≈ ${s.minFeature !== null ? s.minFeature.toPrecision(4) : 'n/a'}, smallest gap ≈ ${s.minGap !== null ? s.minGap.toPrecision(4) : 'n/a'} ${g.units}` : 'the file declares no length unit'));
    rep.assumptions.push('a 2-D layout carries no layer thicknesses or materials: a vertical stack must be supplied for 3-D simulation'); facts.calibrated = !!g.units && !/unitless|px/.test(g.units);
  } else if (g.kind === 'rods') V.push(chk('Rod records', g.rods.length > 0, `${g.rods.length} rods on ${g.area.join(' × ')} ${g.units || 'nm'}`));
  return facts;
}
const SUITE_NAME = { 3: 'atomic and crystal structure', 4: 'nanomaterial geometry and morphology', 5: 'metamaterial architecture and unit-cell design', 9: 'fluid dynamics and porous-media transport', 10: 'heat transfer', 12: 'electromagnetic and photonic engineering', 14: 'structural and nanomechanical engineering', 15: 'acoustic and phononic engineering', 16: 'microstructure and defect evolution', 17: 'atomistic simulation', 18: 'multiscale modelling and homogenisation', 19: 'component and system analysis', 21: 'fabrication and process engineering', 23: 'characterisation and metrology', 29: 'reverse engineering and reconstruction' };
function suitability(g, f) {
  const out = [], S = (suiteId, ok, why) => out.push({ suiteId, suite: SUITE_NAME[suiteId], ok: !!ok, why });
  if (!g) return out;
  if (g.kind === 'atoms') { S(3, true, 'atom species and positions' + (g.cell ? ' with lattice vectors' : '') + ' are available'); S(17, true, g.cell ? 'periodic cell and atom list can seed DFT / MD set-ups' : 'finite cluster: suitable for molecular / cluster calculations; periodic codes need a cell'); S(4, true, 'particle size, facets and surface fraction can be derived from the atom positions'); S(23, !!g.cell, g.cell ? 'a powder diffraction pattern can be computed from the cell' : 'no unit cell: Bragg positions cannot be indexed'); S(14, false, 'an atomistic model is not a continuum body: convert to an envelope mesh or use atomistic mechanics'); }
  else if (g.kind === 'mesh') { const closed = f.watertight && f.manifold; S(4, true, 'surface area, volume (if closed) and dimensions can be measured'); S(5, closed || f.solid, closed ? 'closed surface: relative density and unit-cell analysis are possible' + (f.periodic ? ' and the cell is periodic as meshed' : '; periodic faces do not match, so apply periodicity with care') : f.solid ? 'solid-element mesh' : 'open surface: volume fraction is undefined'); S(14, closed || f.solid, closed || f.solid ? 'closed surface or solid mesh can be volume-meshed / analysed' : 'open or non-manifold surface: repair or re-export before structural analysis'); S(9, closed, closed ? 'a closed boundary can define a flow domain (its interior or exterior)' : 'flow domains need a closed boundary'); S(10, closed || f.solid, closed || f.solid ? 'closed body for conduction analysis' : 'open surface'); S(15, closed || f.solid, closed || f.solid ? 'closed body' + (f.periodic ? ' with periodic faces (Bloch analysis possible)' : '') : 'open surface'); S(21, closed, closed ? 'watertight: suitable for additive-manufacturing preparation' : 'slicers need a watertight mesh'); S(19, closed || f.solid, closed || f.solid ? 'usable as a component body' : 'open surface'); S(17, false, 'a continuum surface holds no atoms'); }
  else if (g.kind === 'voxels') { S(16, true, 'phase map on a regular grid'); S(18, f.phases > 1, f.phases > 1 ? 'representative volume with distinct phases for homogenisation' : 'single phase: nothing to homogenise'); S(9, !!f.porePercolates, f.porePercolates ? 'pore phase percolates: permeability and tortuosity can be estimated' : 'pore phase (id 0) does not percolate: no through-flow'); S(10, !!f.solidPercolates || f.phases > 1, 'effective conductivity can be computed on the grid'); S(14, !!f.solidPercolates, f.solidPercolates ? 'load-bearing solid phase is connected' : 'solid phase is not connected across the volume'); S(23, true, 'phase fractions, chord lengths and specific surface can be measured' + (f.calibrated ? '' : ' (in voxel units until calibrated)')); S(29, true, 'measured or generated volume usable for reconstruction'); S(12, false, 'a full-wave optical model needs smooth interfaces or a layout: convert the phase of interest to a mesh first'); }
  else if (g.kind === 'lattice') { if (f.pointCloud) { S(23, true, 'point measurements can be compared with nominal geometry'); S(29, true, 'input for surface reconstruction'); S(14, false, 'no surface or connectivity'); S(5, false, 'no connectivity'); } else if (f.toolpath) { S(21, true, 'toolpath statistics (length, layers, time) are available'); S(14, false, 'a toolpath is not a solid model'); } else { S(5, true, 'node–edge architecture: connectivity, Maxwell number and unit-cell statistics are available'); S(14, f.connected && f.hasSection && !f.wireframe, f.wireframe ? 'CAD wireframe preview, not a structural model' : !f.hasSection ? 'assign a strut cross-section first' : f.connected ? 'connected beam network with sections' : 'network is not connected'); S(18, !!f.periodic && !f.wireframe, f.periodic ? 'periodic cell: homogenisation is well posed' : 'no periodic match on opposite faces'); S(15, !!f.periodic && !f.wireframe, f.periodic ? 'periodic cell for Bloch-wave analysis' : 'not periodic as given'); S(21, f.hasSection && !f.wireframe, f.hasSection ? 'can be thickened to a printable mesh' : 'no strut diameter'); } }
  else if (g.kind === 'layout') { S(12, true, 'planar pattern for metasurface / resonator analysis (add the layer stack)'); S(21, f.calibrated, f.calibrated ? 'lithography pattern with units: feature and gap sizes can be checked against process limits' : 'no units: critical dimensions cannot be assessed'); S(5, true, 'planar unit-cell pattern'); S(14, false, 'a 2-D layout has no thickness: extrude it first'); }
  else if (g.kind === 'rods') { S(4, true, 'rod array statistics'); S(12, true, 'array geometry for optical / field-emission estimates'); }
  return out;
}

// ---------------------------------------------------------------- kernel and native-container paths
const KERNEL_FMT = { step: 'step', iges: 'iges', brep: 'brep', '3dm': '3dm', ifc: 'ifc' };
const describeDeflection = d => !d ? '' : d.type === 'segments' ? `${d.circleSegments} segments per full circle (curved profiles and sweeps)` : `linear deflection ${d.linear}${d.type === 'absolute' ? ' mm' : ' × the mean bounding-box edge'}${d.angular ? `, angular deflection ${d.angular} rad` : ''}`;
const errText = e => (e && e.message ? e.message : String(e));
// Tessellate through a kernel. Returns {g} (a mesh result carrying notes / lost / meta, like the built-in parsers) or {note}
// saying why the kernel path was not taken.
async function kernelImport(format, ctx, extra) {
  const fid = KERNEL_FMT[format.id], kn = K.KERNELS[K.kernelFor(fid)]; let t;
  try { t = await K.tessellate({ format: fid, bytes: ctx.bytes, options: ctx.opts.tessellation || {} }); } catch (e) { return { note: `the ${kn.engine} kernel (${kn.name}) could not be used — ${errText(e)}` }; }
  if (!t.mesh.faces.length) return { note: `the ${kn.engine} kernel read the file but produced no triangles${t.warnings.length ? ` (${t.warnings.join('; ')})` : ''}`, empty: t };
  const st = t.stats, exact = fid !== '3dm' && fid !== 'ifc', notes = [], lost = [], meta = { ...(t.meta || {}) };
  notes.push(`geometry path: ${kn.engine} kernel (${kn.name} ${kn.version}, loaded on demand) — ${st.bodies} ${st.bodies === 1 ? 'body' : 'bodies'} (${st.solids} closed), ${exact || st.faces ? `${st.faces} B-rep faces, ` : ''}${st.triangles} triangles in ${st.seconds.toFixed(2)} s`);
  if (fid === 'step' || fid === 'iges' || fid === 'brep') {
    notes.push(`the mesh is a tessellation of the exact B-rep at ${describeDeflection(st.deflection)}: curved faces are approximated to within that tolerance and the analytic / NURBS surface definitions are not retained`);
    let light = null, declared = ''; if (fid !== 'brep' && ctx.bytes.length <= 3e7) { try { const lr = fid === 'step' ? F.parseSTEP(ctx.bytes) : F.parseIGES(ctx.bytes); light = lr.meta || null; declared = lr.units || ''; } catch (e) { light = null; } }
    if (light) { Object.assign(meta, light); const u = declared; if (u) meta.declaredUnits = u; notes.push(`coordinates are in millimetres: the kernel converts from the unit the file declares${u ? ` (${u})` : ''}`); const fileFaces = light.entityCounts ? (light.entityCounts.ADVANCED_FACE || 0) + (light.entityCounts.FACE_SURFACE || 0) : 0; if (fid === 'step' && fileFaces > st.faces) lost.push(`the file lists ${fileFaces} faces but the kernel transferred ${st.faces}: some faces were not translated`); }
    else if (fid !== 'brep') notes.push('coordinates are in millimetres: the kernel converts from the unit the file declares');
    else notes.push('a BREP file declares no length unit: coordinates are kept as stored (OpenCASCADE applications such as FreeCAD write millimetres)');
    lost.push(fid === 'brep' ? 'BREP stores shapes only: there are no names, colours or units to import' : 'the product structure is flattened to a list of bodies (name and colour kept per body, instances placed); PMI, annotations, layers and validation properties are not imported');
  } else if (fid === '3dm') { notes.push('3DM geometry comes from what the file stores: mesh objects and saved render meshes are used as they are, and are not re-tessellated to a deflection'); lost.push('materials, textures, groups, views, annotations and curves are not imported'); }
  else { notes.push(`IFC body geometry is tessellated in metres with ${describeDeflection(st.deflection)}; planar faces are exact`); lost.push('property sets, quantities, materials, the spatial structure and element relationships are not imported'); }
  for (const w of t.warnings) (/not tessellated|could not|missing|not converted|no stored|only some/.test(w) ? lost : notes).push(w);
  const g = { ...t.mesh, notes, lost, meta }; if (t.parts.length > 1) g.shells = `${t.parts.length} bodies kept as separate shells (touching bodies are not merged)`;
  if (st.untessellatedFaces || (fid === '3dm' && t.meta && t.meta.converted && (t.meta.converted.brepNone || t.meta.converted.brepPartial || t.meta.converted.extrusionNone))) g.complete = false;
  extra.parts = t.parts; extra.kernel = { name: kn.name, engine: kn.engine, version: kn.version, licence: kn.licence, stats: st };
  return { g };
}
// FreeCAD: every stored shape is an OpenCASCADE BREP, tessellated one object at a time and placed by its container.
async function fcstdImport(ctx, extra) {
  const d = N.readFCStd(ctx.bytes), kn = K.KERNELS.occt; extra.meta = { ...(extra.meta || {}), ...d.meta }; if (d.extra.preview) extra.preview = d.extra.preview;
  const info = [`FreeCAD document${d.meta.programVersion ? ` written by FreeCAD ${d.meta.programVersion}` : ''}: ${d.meta.objects.length} objects, ${d.shapes.length} with a stored shape to import`]; if (extra.preview) info.push(`thumbnail (${extra.preview.mime}, ${extra.preview.bytes.length} bytes)`);
  if (!d.shapes.length) return { info, note: 'the document stores no visible solid or surface shapes', notes: d.notes };
  const vertices = [], faces = [], parts = [], notes = [...d.notes], lost = []; let failed = 0, secs = 0, brepFaces = 0, defl = null;
  for (const s of d.shapes) {
    let t; try { t = await K.tessellate({ format: 'brep', bytes: s.brep, options: ctx.opts.tessellation || {} }); } catch (e) { if (!parts.length && !failed && /could not be loaded|no kernel/.test(errText(e))) return { info, note: `the ${kn.engine} kernel (${kn.name}) could not be used — ${errText(e)}`, notes: d.notes }; failed++; lost.push(`object "${s.label}" (${s.type}): the kernel could not read its shape (${errText(e)})`); continue; }
    secs += t.stats.seconds; defl = t.stats.deflection; if (!t.mesh.faces.length) continue; const base = vertices.length, f0 = faces.length, M = s.transform;
    for (const p of t.mesh.vertices) vertices.push(M ? [M[0] * p[0] + M[1] * p[1] + M[2] * p[2] + M[3], M[4] * p[0] + M[5] * p[1] + M[6] * p[2] + M[7], M[8] * p[0] + M[9] * p[1] + M[10] * p[2] + M[11]] : p);
    for (const f of t.mesh.faces) faces.push([f[0] + base, f[1] + base, f[2] + base]);
    for (const p of t.parts) { parts.push({ ...p, name: t.parts.length > 1 ? `${s.label} / ${p.name}` : s.label, colour: s.colour || p.colour, type: s.type, source: s.file, faceRange: [p.faceRange[0] + f0, p.faceRange[1] + f0] }); brepFaces += p.brepFaces || 0; }
    for (const w of t.warnings) (/not tessellated|could not|missing/.test(w) ? lost : notes).push(`${s.label}: ${w}`);
  }
  if (!faces.length) return { info, note: 'none of the stored shapes produced triangles (wires, sketches or empty shapes only)', notes: d.notes };
  notes.unshift(`geometry path: ${kn.engine} kernel (${kn.name} ${kn.version}, loaded on demand) — ${parts.length} ${parts.length === 1 ? 'body' : 'bodies'} (${parts.filter(p => p.closed).length} closed) from ${d.shapes.length - failed} FreeCAD ${d.shapes.length - failed === 1 ? 'object' : 'objects'}, ${brepFaces} B-rep faces, ${faces.length} triangles in ${secs.toFixed(2)} s`, `the mesh is a tessellation of the stored OpenCASCADE shapes at ${describeDeflection(defl)}: the exact surfaces are not retained`, 'FreeCAD stores shapes in millimetres');
  lost.push('the parametric history (sketches, constraints, feature parameters, spreadsheets) is not imported: only the resulting shapes are');
  const g = { kind: 'mesh', vertices, faces, units: 'mm', notes, lost, meta: extra.meta }; if (parts.length > 1) g.shells = `${parts.length} bodies kept as separate shells (touching bodies are not merged)`; if (failed) g.complete = false;
  extra.parts = parts; extra.kernel = { name: kn.name, engine: kn.engine, version: kn.version, licence: kn.licence, stats: { bodies: parts.length, solids: parts.filter(p => p.closed).length, faces: brepFaces, triangles: faces.length, seconds: secs, deflection: defl } };
  return { g, info };
}
// Parasolid: each transmit stream is read node by node, translated to STEP (xt.js) and tessellated by OpenCASCADE.
// streams: [{name, bytes}]. Returns {g, info} or {note} when nothing could be translated (the caller falls back to the inventory).
const PS_LOST = 'feature history, sketches, PMI and assembly constraints are not recovered: a Parasolid transmit holds only the resulting boundary representation';
async function parasolidImport(streams, ctx, extra, what) {
  const kn = K.KERNELS.occt, vertices = [], faces = [], parts = [], notes = [], lost = [], info = [], sum = { streams: 0, bodies: 0, faces: 0, translated: 0, edges: 0, approximatedEdges: 0, skipped: {}, approximated: {}, untessellated: 0, schemas: [] }, steps = []; let secs = 0, defl = null, firstError = '';
  for (const st of streams) {
    let model, tr;
    try { model = XT.readXT(st.bytes); tr = XT.xtToStep(model, { name: st.label || st.name || ctx.name }); } catch (e) { if (!firstError) firstError = errText(e); lost.push(`${st.name || 'Parasolid stream'}: the node records could not be translated (${errText(e)})`); continue; }
    const s = tr.stats; sum.streams++; sum.approximatedFaces = (sum.approximatedFaces || 0) + (s.approximatedFaces || 0); sum.bodies += s.bodies; sum.faces += s.faces; sum.translated += s.translated; sum.edges += s.edges; sum.approximatedEdges += s.approximatedEdges; if (!sum.schemas.includes(model.schema)) sum.schemas.push(model.schema);
    for (const [k, v] of Object.entries(s.skipped)) sum.skipped[k] = (sum.skipped[k] || 0) + v; for (const [k, v] of Object.entries(s.approximated)) sum.approximated[k] = (sum.approximated[k] || 0) + v;
    info.push(`Parasolid ${model.encoding} transmit${model.partition ? ' (partition)' : ''}${streams.length > 1 ? ` "${st.name}"` : ''}, modeller version ${model.modellerVersion}, schema ${model.schema}: ${model.nodes.size} node records read (${Object.entries(model.counts).filter(([k]) => ['BODY', 'FACE', 'EDGE', 'VERTEX'].includes(k)).map(([k, v]) => `${v} ${k.toLowerCase()}`).join(', ')})`);
    for (const w of model.warnings) notes.push(w); lost.push(...tr.lost); if (s.wireBodies) lost.push(`${s.wireBodies} wire ${s.wireBodies === 1 ? 'body' : 'bodies'} (curves without faces) ${s.wireBodies === 1 ? 'was' : 'were'} not translated`);
    if (!tr.step) continue; const bytes = new TextEncoder().encode(tr.step); let t;
    try { t = await K.tessellate({ format: 'step', bytes, options: ctx.opts.tessellation || {} }); } catch (e) { if (!parts.length && /could not be loaded|no kernel/.test(errText(e))) return { info, note: `the Parasolid B-rep was translated to STEP (${s.translated} of ${s.faces} faces), but the ${kn.engine} kernel (${kn.name}) could not be used to tessellate it — ${errText(e)}` }; lost.push(`${st.name || 'Parasolid stream'}: the kernel could not read the translated STEP (${errText(e)})`); continue; }
    steps.push({ name: `${String(st.label || st.name || ctx.name || 'parasolid').replace(/\.[^.]+$/, '')}.step`, mime: 'model/step', format: 'step', bytes, description: `STEP AP214 translation of the Parasolid B-rep (${s.translated} of ${s.faces} faces), as handed to the kernel` });
    secs += t.stats.seconds; defl = t.stats.deflection; sum.untessellated += t.stats.untessellatedFaces || 0; if (!t.mesh.faces.length) continue;
    const base = vertices.length, f0 = faces.length, written = tr.bodies.filter(b => b.written); for (const p of t.mesh.vertices) vertices.push(p); for (const f of t.mesh.faces) faces.push([f[0] + base, f[1] + base, f[2] + base]);
    t.parts.forEach((p, i) => { const b = t.parts.length === written.length ? written[i] : written.find(x => x.name === p.name); parts.push({ ...p, colour: p.colour || (b && b.colour) || null, type: b ? `Parasolid ${b.type} body` : 'Parasolid body', source: st.name || '', faceRange: [p.faceRange[0] + f0, p.faceRange[1] + f0] }); });
    for (const w of t.warnings) (/not tessellated|could not|missing/.test(w) ? lost : notes).push(w);
  }
  if (!faces.length) return { info, note: sum.streams ? `the Parasolid data was read (${sum.faces} faces in ${sum.bodies} bodies) but no face could be translated and tessellated${Object.keys(sum.skipped).length ? ` (${Object.entries(sum.skipped).map(([k, v]) => `${v} × ${k}`).join(', ')})` : ''}` : `the Parasolid node records could not be translated (${firstError || 'no transmit stream'})`, lost };
  const closed = parts.filter(p => p.closed).length, nskip = sum.faces - sum.translated;
  notes.unshift(`geometry path: Parasolid B-rep translated to STEP and tessellated by OpenCASCADE (${kn.name} ${kn.version}, loaded on demand) — ${what}; ${parts.length} ${parts.length === 1 ? 'body' : 'bodies'} (${closed} closed), ${sum.translated} of ${sum.faces} faces translated, ${faces.length} triangles in ${secs.toFixed(2)} s`,
    `the mesh is a tessellation of the translated B-rep at ${describeDeflection(defl)}: the analytic and NURBS surface definitions are written to STEP exactly and are not retained in the mesh`,
    'Parasolid stores lengths in metres; the mesh is in millimetres (converted by the kernel from the metre unit declared in the STEP text)');
  if (sum.approximatedEdges) notes.push(`${sum.approximatedEdges} of ${sum.edges} edges lie on curves that have no STEP form (${Object.entries(sum.approximated).map(([k, v]) => `${v} × ${k}`).join(', ')}): each was written as a cubic B-spline through the curve's stored chart points or sampled points, an approximation of the exact curve and not the curve itself`);
  if (sum.approximatedFaces) notes.push(`${sum.approximatedFaces} of ${sum.faces} faces lie on rolling-ball blend surfaces, which have no STEP form: each was written as a bicubic B-spline surface through circular sections of the blend (exact contact points on the supporting surfaces), an approximation of the blend and not the surface itself`);
  if (nskip) lost.push(`${nskip} of ${sum.faces} faces were not translated and are missing from the mesh (${Object.entries(sum.skipped).map(([k, v]) => `${v} × ${k}`).join(', ')}): bodies that lose a face are written as open shells and report no volume`);
  lost.push(PS_LOST);
  const g = { kind: 'mesh', vertices, faces, units: 'mm', notes, lost, meta: {} }; if (parts.length > 1) g.shells = `${parts.length} bodies kept as separate shells (touching bodies are not merged)`; if (nskip || sum.untessellated) g.complete = false;
  extra.parts = parts; extra.parasolid = { ...sum, faces: sum.faces, translated: sum.translated, closedBodies: closed, units: 'm' }; extra.embedded = [...(extra.embedded || []), ...steps];
  extra.kernel = { name: kn.name, engine: kn.engine, version: kn.version, licence: kn.licence, stats: { bodies: parts.length, solids: closed, faces: sum.translated, untessellatedFaces: sum.untessellated, triangles: faces.length, seconds: secs, deflection: defl } };
  return { g, info };
}
// Container-level readers for native CAD files: what they recover is reported item by item.
const NATIVE = { solidworks: c => N.parseSolidWorks(c.bytes), inventor: c => N.parseOleCad(c.bytes, 'Inventor'), solidedge: c => N.parseOleCad(c.bytes, 'Solid Edge'), revit: c => N.parseOleCad(c.bytes, 'Revit'), fusion: c => N.parseFusion(c.bytes), jt: c => N.parseJT(c.bytes), parasolid: c => N.parseParasolid(c.bytes), sat: c => N.parseSAT(c.bytes), nx: c => N.parseNX(c.bytes) };

// ---------------------------------------------------------------- import
// importFile({name, bytes}, opts) → {format, detection, geometry|null, extra, report}.
// opts: format (force an id), repair (default true), files ({name: bytes} companions for .xdmf / .pvd), raw (companion
// data for detached NRRD / MHD / NIfTI headers), nx, ny, nz, dtype (headerless volumes), threshold, dataset (HDF5 path),
// kernel (false: never load a geometry kernel), tessellation ({linearDeflection, angularDeflection, deflectionType}).
// extra.parts lists the bodies of a kernel import ({name, colour, faceRange, area, volume, …}), extra.preview an embedded
// preview image ({mime, bytes}), extra.embedded any neutral data streams unpacked from a native container.
export async function importFile({ name, bytes }, opts = {}) {
  let b = F.toBytes(bytes ?? ''), det = opts.format ? { id: opts.format, confidence: 1, why: 'format chosen by the user' } : detectFormat(name, b);
  const report = { imported: [], assumptions: [], lost: [], validation: [], repairs: [], suitability: [], provenance: { file: name || '', bytes: b.length, format: '', version: null, importedAt: null, units: null, detection: det.why, confidence: det.confidence, conversions: [] } }, extra = {};
  for (let depth = 0; depth < 3 && ['gzip', 'compress', 'zipwrap'].includes(det.id); depth++) {
    const kind = det.id;
    try {
      let innerName = String(name || '').replace(/\.(gz|z)$/i, '').replace(/\.stpz$/i, '.stp');
      if (kind === 'gzip') b = F.gunzip(b); else if (kind === 'compress') b = N.unLZW(b); else { const en = F.unzip(b, WRAPPED).filter(e => e.data && WRAPPED.test(e.name)).sort((x, y) => y.size - x.size)[0]; if (!en) throw new Error('the wrapped file could not be inflated'); b = en.data; innerName = en.name; }
      const inner = detectFormat(innerName, b), label = kind === 'gzip' ? 'gzip container' : kind === 'compress' ? 'Unix compress (.Z) container' : `zip archive (entry ${innerName})`;
      report.imported.push(`${label} inflated to ${b.length} bytes`); report.provenance.conversions.push(kind === 'gzip' ? 'gzip decompression' : kind === 'compress' ? 'LZW (.Z) decompression' : 'zip extraction'); det = { ...inner, why: `${label}; inner content: ${inner.why}` };
    } catch (e) { det = { id: 'unknown', confidence: 0.2, why: `${kind === 'gzip' ? 'gzip' : kind === 'compress' ? 'Unix compress (.Z)' : 'zip'} container that could not be unpacked (${errText(e)})` }; }
    report.provenance.detection = det.why;
  }
  const format = BY[det.id] || UNKNOWN, ctx = { name: name || '', bytes: b, opts };
  report.provenance.format = format.name; report.imported.push(`recognised as ${format.name} (${det.why}; confidence ${(det.confidence * 100).toFixed(0)} %)`);
  if (det.confidence < 0.6 && format.id !== 'unknown') report.assumptions.push('format identification rests on weak evidence (mainly the file extension): confirm that the result is what you expect');
  let g = null, failed = null, native = false, kernelNote = '';
  // 1. exact-geometry kernels (on demand)
  if (opts.kernel !== false && (KERNEL_FMT[format.id] || format.id === 'fcstd')) {
    try {
      if (format.id === 'fcstd') { const r = await fcstdImport(ctx, extra); report.imported.push(...(r.info || [])); if (r.g) g = r.g; else { kernelNote = r.note; report.assumptions.push(...(r.notes || [])); native = true; } }
      else { const r = await kernelImport(format, ctx, extra); if (r.g) g = r.g; else kernelNote = r.note; }
    } catch (e) { kernelNote = `the kernel path failed (${errText(e)})`; }
  } else if (opts.kernel === false && (KERNEL_FMT[format.id] || format.id === 'fcstd')) kernelNote = 'the geometry kernel was switched off for this import';
  // 2. container-level readers for native CAD files
  if (!g && NATIVE[format.id]) {
    try { const r = NATIVE[format.id](ctx); native = true; extra.meta = { ...(extra.meta || {}), ...r.meta }; Object.assign(extra, r.extra); report.imported.push(...r.recovered.map(s => `recovered — ${s}`)); report.lost.push(...r.notRecovered); if (r.geometry) g = r.geometry; }
    catch (e) { report.assumptions.push(`the container could not be opened in detail (${errText(e)}): only the signature-level facts are reported`); }
  }
  // 2b. Parasolid B-reps — a transmit file, or the partitions unpacked from a native container: XT → STEP → OpenCASCADE
  {
    const all = format.id === 'parasolid' ? [{ name: ctx.name, bytes: ctx.bytes }] : (extra.embedded || []).filter(e => e.format === 'parasolid'), ps = format.id === 'solidworks' ? all.slice(0, 1) : all;
    if (ps.length && opts.kernel === false) kernelNote = 'the geometry kernel was switched off for this import';
    else if (ps.length) {
      const what = format.id === 'parasolid' ? 'transmit file' : format.id === 'solidworks' ? `partition of ${/configuration (\d+)/.test(ps[0].description || '') ? /(configuration \d+(?: "[^"]*")?)/.exec(ps[0].description)[1] : 'the first configuration'} unpacked from the SOLIDWORKS container` : `${ps.length} partition stream(s) located in the ${format.name} file`;
      try {
        const r = await parasolidImport(ps.map(e => ({ ...e, label: ctx.name })), ctx, extra, what); report.imported.push(...(r.info || []));
        if (r.g) { g = r.g; native = false; if (all.length > ps.length) report.assumptions.push(`the file holds ${all.length} configuration partitions: the first was tessellated, the others remain available under extra.embedded`); }
        else { kernelNote = r.note; report.lost.push(...(r.lost || [])); }
      } catch (e) { kernelNote = `the Parasolid translation failed (${errText(e)})`; }
    }
  }
  // 3. built-in readers
  if (!g && H[format.id]) { try { g = H[format.id](ctx); } catch (e) { failed = e && e.message ? e.message : String(e); if (e instanceof TypeError || e instanceof RangeError) failed = `${format.name}: the file is truncated or malformed (the reader stopped with: ${failed})`; } }
  if (kernelNote) { const got = !!(g && g.kind), fallback = got ? 'the built-in lightweight reader was used instead' : H[format.id] ? 'the built-in lightweight reader was tried instead' : 'only the inventory below is available'; (got ? report.assumptions : report.lost).push(`geometry path: ${kernelNote}; ${fallback}`); }
  else if (g && extra.kernel) report.provenance.conversions.push(`B-rep tessellation by ${extra.kernel.name} ${extra.kernel.version}`);
  if (failed) { report.validation.push(chk('File could be parsed', false, failed)); report.lost.push('no geometry imported because the file could not be parsed'); if (format.neutral) report.advice = format.neutral; extra.meta = { ...(extra.meta || {}), ...sniff(format.id, ctx) }; }
  if (g && g.kind === undefined) g = null;
  const res = g; if (g && g.kind === null) g = null;
  if (res) {
    report.assumptions.push(...(res.notes || [])); report.lost.push(...(res.lost || []));
    for (const k of ['truncated', 'capped']) if (res[k]) report.lost.push(res[k]);
    for (const k of ['meta', 'frames', 'regions', 'layers', 'cellData', 'orientation', 'quality', 'phases', 'layerStack', 'wireframe', 'sets', 'stats', 'groups', 'beams', 'wire', 'elementCounts']) if (res[k] !== undefined && !(Array.isArray(res[k]) && !res[k].length)) extra[k === 'stats' ? 'fileStats' : k] = k === 'meta' ? { ...(extra.meta || {}), ...res.meta } : res[k];
    if (res.grey) extra.fields = { grey: res.grey, range: res.valueRange, threshold: res.threshold };
    if (res.layerStack) report.imported.push(`layer stack of ${res.layerStack.length} layers`);
  }
  if (!g && !failed) {
    if (!native || format.id === 'fcstd') extra.meta = { ...sniff(format.id, ctx), ...(extra.meta || {}) };
    const facts = native && format.id !== 'fcstd' ? [] : Object.entries(extra.meta).filter(([k, v]) => v !== null && v !== undefined && typeof v !== 'object' && k !== 'bytes').slice(0, 8).map(([k, v]) => `${k}: ${v}`);
    if (facts.length) report.imported.push(`header metadata — ${facts.join('; ')}`);
    if (!res || !res.layerStack) { report.lost.push(format.id === 'unknown' ? 'nothing imported' : native && format.id !== 'fcstd' ? `no geometry reconstructed — ${format.note}` : kernelNote ? `geometry not reconstructed: ${format.name} needs the geometry kernel for its surfaces, and the kernel was not used for this import` : `geometry not reconstructed: ${format.note}`); report.advice = format.neutral || UNKNOWN.neutral; report.validation.push(chk('Geometry reconstructed', false, format.support === 'inventory' ? `${format.name} is supported at inventory level only (recognition and header metadata)` : native ? `${format.name}: metadata was recovered from the container, but its geometry cannot be rebuilt without the vendor's kernel` : kernelNote ? `${format.name}: the geometry kernel was not available and there is no built-in reader for this geometry` : 'the file was recognised but holds no geometry that this importer can rebuild')); }
  }
  if (g) {
    for (const k of ['notes', 'lost']) delete g[k];
    if (extra.kernel) delete g.meta;
    const n = g.kind === 'mesh' ? `${g.vertices.length} vertices and ${g.faces.length} triangles` : g.kind === 'atoms' ? `${g.atoms.length} atoms${g.cell ? ' in a periodic cell' : ''}` : g.kind === 'voxels' ? `${g.nx} × ${g.ny} × ${g.nz} voxels` : g.kind === 'lattice' ? (g.pointCloud ? `${g.nodes.length} points` : `${g.nodes.length} nodes and ${g.edges.length} edges`) : g.kind === 'layout' ? `${g.polygons.length} polygons` : `${(g.rods || []).length} rods`;
    report.imported.push(`${g.kind} geometry: ${n}`); if (extra.frames) report.imported.push(`${extra.frames.length} time frames for animation`); if (extra.elementCounts) report.imported.push(`elements — ${Object.entries(extra.elementCounts).map(([k, v]) => `${k}: ${v}`).join(', ')}`);
    if (extra.parts) report.imported.push(`bodies — ${extra.parts.slice(0, 12).map(p => `${p.name}${p.closed && p.volume !== undefined ? ` (closed, ${p.volume.toPrecision(5)} ${g.units || 'units'}³)` : ' (open surface)'}`).join('; ')}${extra.parts.length > 12 ? `; … ${extra.parts.length - 12} more` : ''}`);
    if (extra.frames) g.frames = extra.frames;
    const nf = g.kind === 'mesh' ? g.faces.length : 0, facts = validate(g, report, extra, opts); report.suitability = suitability(g, facts);
    if (extra.parts && g.kind === 'mesh' && g.faces.length !== nf) { for (const p of extra.parts) p.faceRange = null; report.assumptions.push('the mesh repair changed the triangle list, so the per-body triangle ranges no longer apply and were cleared'); }
    report.provenance.units = g.units || null; if (!g.units) report.assumptions.push('the file declares no length unit');
    if (g.complete === false) report.validation.push(chk('Complete reconstruction', false, 'only part of the file\'s geometry could be rebuilt (see "lost")'));
    if (native && format.neutral && g.wireframe) report.advice = format.neutral;
  }
  if (extra.embedded) report.assumptions.push(`${extra.embedded.length} embedded data stream(s) are available under extra.embedded for saving`);
  report.provenance.version = (extra.meta && (extra.meta.version ?? extra.meta.protocol ?? extra.meta.schema ?? extra.meta.release ?? extra.meta.igesVersionFlag ?? extra.meta.superblockVersion)) ?? null;
  report.provenance.repairs = report.repairs.length; report.provenance.support = format.support; report.provenance.geometryPath = extra.kernel ? `kernel: ${extra.kernel.name} ${extra.kernel.version}` : g ? 'built-in reader' : 'none';
  report.assumptions = [...new Set(report.assumptions)]; report.lost = [...new Set(report.lost)];
  return { format, detection: det, geometry: g, extra, report };
}

// ---------------------------------------------------------------- representation conversion
const VDW = { H: 1.2, He: 1.4, Li: 1.82, C: 1.7, N: 1.55, O: 1.52, F: 1.47, Ne: 1.54, Na: 2.27, Mg: 1.73, Al: 1.84, Si: 2.1, P: 1.8, S: 1.8, Cl: 1.75, Ar: 1.88, K: 2.75, Ca: 2.31, Ti: 2.15, Fe: 2.04, Ni: 1.63, Cu: 1.4, Zn: 1.39, Ga: 1.87, Ge: 2.11, As: 1.85, Se: 1.9, Br: 1.85, Pd: 1.63, Ag: 1.72, Cd: 1.58, In: 1.93, Sn: 2.17, Te: 2.06, I: 1.98, Pt: 1.75, Au: 1.66, Pb: 2.02 };
function meshToVoxels(m, { resolution = 64 } = {}) { // ray parity along x through voxel centres
  const [lo, hi] = bboxOf(m.vertices), ext = [0, 1, 2].map(a => max(1e-12, hi[a] - lo[a])), res = max(4, min(256, round(resolution))), h = max(...ext) / res, d = ext.map(e => max(1, Math.ceil(e / h)));
  const [nx, ny, nz] = d, data = new Uint8Array(nx * ny * nz), rows = new Array(ny * nz), V = m.vertices, e1 = 1e-7 * h;
  for (const f of m.faces) {
    const A0 = V[f[0]], B = V[f[1]], C = V[f[2]], y0 = min(A0[1], B[1], C[1]), y1 = max(A0[1], B[1], C[1]), z0 = min(A0[2] || 0, B[2] || 0, C[2] || 0), z1 = max(A0[2] || 0, B[2] || 0, C[2] || 0);
    const j0 = max(0, floor((y0 - lo[1]) / h - 0.5)), j1 = min(ny - 1, Math.ceil((y1 - lo[1]) / h - 0.5)), k0 = max(0, floor((z0 - lo[2]) / h - 0.5)), k1 = min(nz - 1, Math.ceil((z1 - lo[2]) / h - 0.5));
    const by = B[1] - A0[1], bz = (B[2] || 0) - (A0[2] || 0), cy = C[1] - A0[1], cz = (C[2] || 0) - (A0[2] || 0), den = by * cz - bz * cy; if (abs(den) < 1e-300) continue;
    for (let k = k0; k <= k1; k++) for (let j = j0; j <= j1; j++) { const py = lo[1] + (j + 0.5) * h + e1 - A0[1], pz = lo[2] + (k + 0.5) * h + 2.3 * e1 - (A0[2] || 0), u = (py * cz - pz * cy) / den, v = (by * pz - bz * py) / den; if (u < 0 || v < 0 || u + v > 1) continue; (rows[j + ny * k] || (rows[j + ny * k] = [])).push(A0[0] + u * (B[0] - A0[0]) + v * (C[0] - A0[0])); }
  }
  let odd = 0;
  for (let r = 0; r < rows.length; r++) { const xs = rows[r]; if (!xs) continue; xs.sort((a, b) => a - b); if (xs.length % 2) odd++; for (let q = 0; q + 1 < xs.length; q += 2) { const i0 = max(0, Math.ceil((xs[q] - lo[0]) / h - 0.5)), i1 = min(nx - 1, floor((xs[q + 1] - lo[0]) / h - 0.5)); for (let i = i0; i <= i1; i++) data[i + nx * r] = 1; } }
  return { vox: { kind: 'voxels', nx, ny, nz, data, voxel: h, units: m.units || '', origin: [...lo] }, odd };
}
function stampSpheres(pts, radii, ids, { voxel, maxDim = 96, units }) {
  let lo = [Infinity, Infinity, Infinity], hi = [-Infinity, -Infinity, -Infinity], rmax = 0; pts.forEach((p, i) => { rmax = max(rmax, radii[i]); for (let a = 0; a < 3; a++) { lo[a] = min(lo[a], p[a]); hi[a] = max(hi[a], p[a]); } });
  lo = lo.map(v => v - rmax); hi = hi.map(v => v + rmax); const ext = hi.map((v, a) => v - lo[a]), h = voxel || max(...ext) / maxDim, d = ext.map(e => max(1, min(400, Math.ceil(e / h)))), [nx, ny, nz] = d, data = new Uint8Array(nx * ny * nz);
  pts.forEach((p, i) => { const r = radii[i], R = Math.ceil(r / h), c = p.map((v, a) => (v - lo[a]) / h - 0.5), ci = c.map(round); for (let dz = -R; dz <= R; dz++) for (let dy = -R; dy <= R; dy++) for (let dx = -R; dx <= R; dx++) { const x = ci[0] + dx, y = ci[1] + dy, z = ci[2] + dz; if (x < 0 || y < 0 || z < 0 || x >= nx || y >= ny || z >= nz) continue; if (((x - c[0]) ** 2 + (y - c[1]) ** 2 + (z - c[2]) ** 2) * h * h <= r * r) data[x + nx * (y + ny * z)] = ids[i]; } });
  return { kind: 'voxels', nx, ny, nz, data, voxel: h, units, origin: lo };
}
const STEP_FN = {
  'atoms>voxels': (g, o) => { const els = [...new Set(g.atoms.map(a => a.el))].sort(), at = g.atoms.slice(0, 200000), v = stampSpheres(at.map(a => [a.x, a.y, a.z]), at.map(a => VDW[a.el] || 2.0), at.map(a => els.indexOf(a.el) + 1), { voxel: o.voxel_A, maxDim: o.resolution || 96, units: 'Å' }); v.phaseLegend = Object.fromEntries(els.map((e, i) => [i + 1, e])); return { g: v, lost: [`atoms replaced by van der Waals spheres (Bondi-type radii; 2.0 Å where no value is tabulated) on a ${v.voxel.toPrecision(3)} Å grid: positions below the voxel size, bonding and periodic images are lost`, ...(g.atoms.length > at.length ? [`only the first ${at.length} atoms were rasterised`] : []), ...(g.cell ? ['the periodic cell is not represented: the grid is the bounding box of the atoms'] : [])] }; },
  'voxels>mesh': (g, o) => { const src = o.phase !== undefined ? { ...g, data: Uint8Array.from(g.data, v => v === o.phase ? 1 : 0), field: undefined } : g, m = marchingCubes(src); return { g: m, lost: [`iso-surface of ${o.phase !== undefined ? 'phase ' + o.phase : 'all non-zero phases'} by marching tetrahedra: the surface is faceted at the voxel scale (staircase smoothing to half-voxel positions)`, 'phase identities other than solid / void are merged', 'the surface is closed at the volume boundary (cut faces are capped)'] }; },
  'mesh>voxels': (g, o) => { const r = meshToVoxels(g, o); return { g: r.vox, lost: [`solid rasterised by ray parity on a ${r.vox.voxel.toPrecision(3)} ${g.units || 'unit'} grid: features thinner than a voxel vanish and surfaces become staircases`, ...(r.odd ? [`${r.odd} grid rows crossed the surface an odd number of times (open or self-intersecting mesh): the fill is unreliable there`] : [])] }; },
  'lattice>mesh': (g, o) => { const [lo, hi] = bboxOf(g.nodes), need = !(g.radius > 0) && !g.radii, src = need ? { ...g, radius: o.radius || 0.01 * hypot(hi[0] - lo[0], hi[1] - lo[1], hi[2] - lo[2]) } : g, m = latticeToMesh(src, { segments: o.segments || 6 }); return { g: m, lost: ['struts become prisms and nodes icosahedra as overlapping closed shells (no Boolean union, no fillets)', ...(need ? [`no strut radius was defined: ${src.radius.toPrecision(3)} ${g.units || 'units'} (1 % of the bounding-box diagonal) was assumed`] : []), ...(m.truncated ? [m.truncated] : [])] }; },
  'lattice>voxels': (g, o) => { const [lo, hi] = bboxOf(g.nodes), ext = max(hi[0] - lo[0], hi[1] - lo[1], hi[2] - lo[2]) || 1, res = o.resolution || 96, r0 = g.radius > 0 ? g.radius : (o.radius || ext / res), h = ext / res, pts = [], rr = []; g.edges.forEach((e, k) => { const a = g.nodes[e[0]], b = g.nodes[e[1]], L = hypot(b[0] - a[0], b[1] - a[1], (b[2] || 0) - (a[2] || 0)), n = max(1, Math.ceil(L / (0.5 * h))), r = g.radii ? g.radii[k] : r0; for (let i = 0; i <= n && pts.length < 3e6; i++) { const t = i / n; pts.push([a[0] + t * (b[0] - a[0]), a[1] + t * (b[1] - a[1]), (a[2] || 0) + t * ((b[2] || 0) - (a[2] || 0))]); rr.push(max(r, 0.5 * h)); } }); const v = stampSpheres(pts, rr, new Uint8Array(pts.length).fill(1), { maxDim: res, units: g.units || '' }); return { g: v, lost: [`struts rasterised as capsules on a ${v.voxel.toPrecision(3)} ${g.units || 'unit'} grid; struts thinner than one voxel are thickened to one voxel`, ...(g.radius > 0 || g.radii ? [] : ['no strut radius was defined: one-voxel struts were drawn'])] }; },
  'rods>mesh': g => { const m = rodsToMesh(g); return { g: m, lost: ['rods become ideal prisms standing on z = 0 (no tapering, tips or substrate)', ...(m.truncated ? [m.truncated] : [])] }; },
  'layout>mesh': (g, o) => { const t = o.thickness > 0 ? o.thickness : 100, parts = []; let skipped = 0; for (const q of g.polygons.slice(0, 20000)) { const tr = triangulate(q.pts); if (!tr.length) { skipped++; continue; } const z0 = o.layerZ && o.layerZ[q.layer] !== undefined ? o.layerZ[q.layer] : 0; parts.push(extrude(q.pts, tr, z0, z0 + t, g.units || 'nm')); } if (!parts.length) throw new Error('layout→mesh: no polygon could be triangulated'); const m = mergeMeshes(parts, g.units || 'nm'); delete m.ranges; m.shells = 'one closed prism per polygon (abutting or overlapping polygons are not merged)'; return { g: m, lost: [`every polygon extruded to a prism of thickness ${t} ${g.units || 'nm'}${o.thickness > 0 ? '' : ' (default: no thickness is stored in a layout)'}; sidewall angle, rounding and layer-specific heights are not modelled`, 'polygons that touch or overlap remain separate shells', ...(skipped ? [`${skipped} self-intersecting or degenerate polygons skipped`] : []), ...(g.polygons.length > 20000 ? ['only the first 20 000 polygons were extruded'] : [])] }; },
  'mesh>lattice': (g, o) => { const n = g.vertices.length, seen = new Set(), edges = []; for (const f of g.faces) for (let e = 0; e < 3; e++) { const a = f[e], b = f[(e + 1) % 3], k = a < b ? a * n + b : b * n + a; if (!seen.has(k)) { seen.add(k); edges.push([a, b]); } } return { g: { kind: 'lattice', nodes: g.vertices.map(p => [p[0], p[1], p[2] || 0]), edges, radius: o.radius || 0, units: g.units || '' }, lost: ['only the edge graph of the triangulation is kept: faces, enclosed volume and surface orientation are discarded', 'triangulation diagonals appear as struts although they are not physical members', ...(o.radius ? [] : ['no strut radius assigned'])] }; },
};
STEP_FN['atoms>mesh'] = (g, o) => { const a = STEP_FN['atoms>voxels'](g, o), b = STEP_FN['voxels>mesh'](a.g, {}); return { g: b.g, lost: [...a.lost, 'envelope = iso-surface of the van der Waals sphere union: element identity is lost'] }; };
const KINDS = ['atoms', 'lattice', 'mesh', 'voxels', 'layout', 'rods'];
// convert(geometry, targetKind, opts) → {geometry, lost, path}. Multi-step routes (for example rods → mesh → voxels)
// are found automatically and every step's loss is reported.
export function convert(geometry, targetKind, opts = {}) {
  if (!geometry || !KINDS.includes(geometry.kind)) throw new Error(`convert: unknown source geometry kind "${geometry && geometry.kind}"`);
  if (!KINDS.includes(targetKind)) throw new Error(`convert: unknown target kind "${targetKind}" (use one of ${KINDS.join(', ')})`);
  if (geometry.kind === targetKind) return { geometry, lost: [], path: [targetKind] };
  const prev = { [geometry.kind]: null }, queue = [geometry.kind];
  while (queue.length && !(targetKind in prev)) { const k = queue.shift(); for (const key of Object.keys(STEP_FN)) { const [a, b] = key.split('>'); if (a === k && !(b in prev)) { prev[b] = k; queue.push(b); } } }
  if (!(targetKind in prev)) throw new Error(`convert: there is no physically meaningful route from ${geometry.kind} to ${targetKind} (${targetKind === 'atoms' ? 'a continuum shape does not define atom positions: build the crystal in the atomic engine and carve it to this envelope instead' : targetKind === 'layout' ? 'a 3-D shape has no unique 2-D mask: define the section plane and layer explicitly' : targetKind === 'rods' ? 'rod arrays are generated, not derived' : 'not implemented'})`);
  const path = [targetKind]; while (prev[path[0]] !== null) path.unshift(prev[path[0]]);
  let g = geometry; const lost = [];
  for (let i = 0; i + 1 < path.length; i++) { const r = STEP_FN[path[i] + '>' + path[i + 1]](g, opts); g = r.g; lost.push(...r.lost.map(s => (path.length > 2 ? `[${path[i]} → ${path[i + 1]}] ` : '') + s)); }
  return { geometry: g, lost, path };
}

// ---------------------------------------------------------------- export
const E = (name, ext, mime, accepts, write, loses = []) => ({ name, ext, mime, accepts, write, loses });
const T = f => (g, o) => ({ text: f(g, o) }), B = f => (g, o) => ({ bytes: f(g, o) });
export const EXPORTERS = {
  stl: E('STL (binary)', '.stl', 'model/stl', ['mesh'], B((g, o) => F.toSTL(g, { ...o, binary: true })), ['STL stores no units, materials or regions; coordinates are written as 32-bit floats']),
  'stl-ascii': E('STL (ASCII)', '.stl', 'model/stl', ['mesh'], T((g, o) => F.toSTL(g, { ...o, binary: false })), ['STL stores no units, materials or regions']),
  obj: E('Wavefront OBJ', '.obj', 'model/obj', ['mesh'], T(F.toOBJ), ['units only as a comment']), ply: E('PLY', '.ply', 'application/octet-stream', ['mesh'], (g, o = {}) => o.binary ? { bytes: F.toPLY(g, o) } : { text: F.toPLY(g, o) }, ['units only as a comment']),
  off: E('OFF', '.off', 'text/plain', ['mesh'], T(F.toOFF), ['no units']), '3mf': E('3MF', '.3mf', 'model/3mf', ['mesh'], B(F.to3MF), ['single object without materials or colour']), amf: E('AMF', '.amf', 'application/x-amf', ['mesh'], T(F.toAMF), ['single volume without materials']),
  step: E('STEP AP214 (faceted B-rep)', '.step', 'model/step', ['mesh'], T(F.toSTEP), ['the solid is a faceted B-rep of planar triangles: no analytic or NURBS surfaces are created']),
  vtk: E('VTK legacy (ASCII)', '.vtk', 'text/plain', ['mesh', 'voxels', 'lattice'], T(F.toVTK), ['no units']), vtu: E('VTK unstructured grid', '.vtu', 'application/xml', ['mesh', 'lattice'], T(F.toVTU), ['no units']), vti: E('VTK image data', '.vti', 'application/xml', ['voxels'], T(F.toVTI), ['no units']),
  msh: E('Gmsh MSH 2.2', '.msh', 'text/plain', ['mesh', 'lattice'], T(F.toMSH), ['surface triangles or line elements only (no volume mesh is generated)']), inp: E('Abaqus input', '.inp', 'text/plain', ['mesh', 'lattice'], T(F.toINP), ['nodes and S3 / B31 elements only: no sections, materials or steps']), bdf: E('NASTRAN bulk data', '.bdf', 'text/plain', ['mesh', 'lattice'], T(F.toBDF), ['GRID and CTRIA3 / CROD cards only: property and material cards must be added']),
  gds: E('GDSII', '.gds', 'application/octet-stream', ['layout'], B(F.toGDS), ['coordinates snapped to the database grid (default 1 nm); layer thickness and materials are not part of GDSII']), oasis: E('OASIS', '.oas', 'application/octet-stream', ['layout'], B(F.toOASIS), ['coordinates snapped to the grid (default 1 nm); no properties or layer names']),
  dxf: E('DXF', '.dxf', 'image/vnd.dxf', ['layout', 'lattice'], T(F.toDXF), ['minimal DXF without tables or blocks']), svg: E('SVG', '.svg', 'image/svg+xml', ['layout'], T(F.toSVG), ['display colours are arbitrary; units recorded in a data attribute only']),
  hdf5: E('HDF5', '.h5', 'application/x-hdf5', ['voxels'], B(F.toHDF5), ['phase ids only; auxiliary fields are not written']), tiff: E('TIFF stack', '.tif', 'image/tiff', ['voxels'], B(F.toTIFFStack), ['phase ids stored as grey values; calibration follows the ImageJ convention']), raw: E('RAW volume', '.raw', 'application/octet-stream', ['voxels'], B(F.toRAW), ['no header at all: note the dimensions, uint8 type and x-fastest order separately']),
  nrrd: E('NRRD', '.nrrd', 'application/octet-stream', ['voxels'], B(F.toNRRD), []), mha: E('MetaImage', '.mha', 'application/octet-stream', ['voxels'], B(F.toMHA), ['units not stored (MetaImage convention is millimetres)']), xdmf: E('XDMF (inline data)', '.xdmf', 'application/xml', ['mesh', 'voxels', 'lattice'], T(F.toXDMF), ['no units']),
  scad: E('OpenSCAD', '.scad', 'text/plain', ['mesh', 'lattice', 'rods', 'layout'], T(F.toSCAD), ['explicit geometry, not a parametric script']), geo: E('Gmsh geometry', '.geo', 'text/plain', ['lattice', 'layout', 'mesh'], T(F.toGEO), ['explicit points and lines without parameters']), json: E('JSON (native)', '.json', 'application/json', KINDS, T(F.toJSON), []),
  cif: E('CIF', '.cif', 'chemical/x-cif', ['atoms'], T(A.toCIF), []), xyz: E('XYZ', '.xyz', 'chemical/x-xyz', ['atoms'], T(A.toXYZ), ['no cell or periodicity']), extxyz: E('Extended XYZ', '.extxyz', 'chemical/x-xyz', ['atoms'], T(A.toExtXYZ), []), poscar: E('VASP POSCAR', '.vasp', 'text/plain', ['atoms'], T(A.toPOSCAR), ['needs a periodic cell']),
  xsf: E('XSF', '.xsf', 'text/plain', ['atoms'], T(A.toXSF), []), pdb: E('PDB', '.pdb', 'chemical/x-pdb', ['atoms'], T(A.toPDB), ['fixed-column precision of 0.001 Å']), lammps: E('LAMMPS data', '.data', 'text/plain', ['atoms'], T(A.toLAMMPSData), ['atom types numbered by element; no force-field terms']), mol: E('MDL MOL', '.mol', 'chemical/x-mdl-molfile', ['atoms'], T(A.toMOL), ['V2000 limit of 999 atoms; no cell']),
};
// exportGeometry(geometry, formatId, opts) → {name, mime, text|bytes, lost}. Converts first when the format cannot hold
// the geometry's representation (for example voxels → STL goes through an iso-surface) and reports what that costs.
export function exportGeometry(geometry, formatId, opts = {}) {
  const ex = EXPORTERS[formatId]; if (!ex) throw new Error(`exportGeometry: unknown format "${formatId}". Available: ${Object.keys(EXPORTERS).join(', ')}`);
  if (!geometry || !geometry.kind) throw new Error('exportGeometry: no geometry to export');
  let g = geometry; const lost = [];
  if (!ex.accepts.includes(g.kind)) { let done = null, err = ''; for (const k of ex.accepts) { try { done = convert(g, k, opts); break; } catch (e) { err = e.message; } } if (!done) throw new Error(`exportGeometry: ${ex.name} holds ${ex.accepts.join(' / ')} geometry and a ${g.kind} geometry cannot be converted to it (${err})`); g = done.geometry; lost.push(`converted ${done.path.join(' → ')} before writing`, ...done.lost); }
  const out = ex.write(g, opts);
  return { name: (opts.name || 'geometry').replace(/\.[A-Za-z0-9]+$/, '') + ex.ext, mime: ex.mime, ...out, lost: [...lost, ...ex.loses] };
}
