// Creation environments and conditions knowledge base.
// Typical windows for the 146 creation environments (specification categories A–I), the condition categories
// (areas 30–34), the metamaterial-class routes (table J) and an explainable rule-based recommender.
// Windows are typical engineering ranges from standard process practice, not limits of any specific tool.

const A = [9e4, 1.1e5];            // ambient pressure window, Pa
const L = s => (s ? s.split('|') : []);
export const ENV_CONDITIONS = {};
export const ENV_TAGS = {};        // compact capability tags used by recommendEnvironments (see TAG_KEY)

// Tag key: N nanomaterial, M metamaterial; 0/1/2/3 product dimensionality (first = principal form);
// me ox sc ca po ce co bi material families (first = principal); al alignment, cf conformal, pe periodic,
// lc low-cost, 3a 3-D architecture, hp high-purity, ep epitaxial; L/P/I lab, pilot, industrial scale;
// F enabling facility or condition combined with other routes; X natural formation environment.
export const TAG_KEY = { N: 'nanomaterial', M: 'metamaterial', me: 'metal', ox: 'oxide', sc: 'semiconductor', ca: 'carbon', po: 'polymer', ce: 'ceramic', co: 'composite', bi: 'biological', al: 'alignment', cf: 'conformal', pe: 'periodic', lc: 'low-cost', '3a': '3D-architecture', hp: 'high-purity', ep: 'epitaxial', L: 'lab', P: 'pilot', I: 'industrial', F: 'enabling facility', X: 'natural' };

function e(id, tags, T_C, P_Pa, time, d, chemical, thermal, transport, fields, mechanical, controls, outputsControlled, suites, hazards, sustainability, resolution_nm, rate) {
  ENV_TAGS[id] = tags;
  ENV_CONDITIONS[id] = {
    T_C, P_Pa, time, atmosphere: d[0], medium: d[1], energy: d[2], scale: d[3],
    chemical: L(chemical), thermal: L(thermal), transport: L(transport), fields: L(fields), mechanical: L(mechanical),
    controls: L(controls), outputsControlled: L(outputsControlled), suites, hazards: L(hazards), sustainability,
    resolution_nm: resolution_nm ?? null, rate: rate ?? null
  };
}

// ── A. Chemical and solution-based environments ──
e(1,'N 1 0 3 ox sc ce al lc L P I',[80,250],[1e5,4e6],[3600,172800],
['sealed aqueous solution, autogenous pressure','liquid water','resistive oven or band heater; optional electric-field or microwave assist','mL lab autoclaves to m³ industrial vessels'],
'precursor 1–100 mM (e.g. Zn(NO₃)₂ + HMTA)|pH 6–12|supersaturation set by slow base release','heating rate 1–10 K/min|isothermal hold 1–48 h|autogenous pressure = water saturation pressure','diffusion-limited in unstirred vessels|fill fraction 50–80 %|precursor depletion limits length',
'optional DC field 0.1–10 kV/m','vessel rated above saturation pressure|seeded substrate facing down avoids sediment',
'temperature|time|precursor concentration|pH','length|diameter|aspect ratio|crystallinity',[6,7,9,10,21,27],'pressurised hot vessel|caustic or amine solutions',
'Water solvent and low temperature; main burdens are heating energy and nitrate/amine effluent — recover heat and treat effluent.',null,'growth 1–10 nm/min');
e(2,'N 0 1 2 sc ox me hp L P',[100,350],[1e5,2e7],[3600,259200],
['sealed organic solvent, autogenous pressure','non-aqueous liquid (ethanol, ethylene glycol, DMF, oleylamine)','oven or heating mantle','mL autoclaves to 100 L-class vessels'],
'precursor 1–200 mM|solvent polarity and boiling point set solubility|capping ligands (oleic acid, PVP)','heating rate 1–10 K/min|hold 2–72 h|pressure from solvent vapour and decomposition gases','natural convection or stirring|fill fraction 40–70 %',
'','liner compatible with solvent (PTFE ≤ 220 °C, PPL ≤ 280 °C)|burst disc sized for gas evolution',
'solvent|temperature|ligand type and ratio|time','size|phase|shape|size distribution',[6,7,8,10,17,21],'flammable solvent under pressure|toxic decomposition gases',
'Organic solvents dominate the footprint; choose recoverable low-toxicity solvents (ethanol, glycols) and distil for reuse.',null,'batch yield g–kg per run');
e(3,'N 0 3 2 ox ce co cf lc L P I',[20,90],A,[3600,604800],
['air, controlled humidity','alcohol–water sol turning to gel','ambient or mild heating; calcination afterwards','beaker to tonne-scale coating and powder plants'],
'alkoxide 0.1–2 M (TEOS, Ti-isopropoxide)|water:alkoxide ratio 1–20|acid or base catalyst pH 1–11','gelation 20–80 °C|drying 60–120 °C|calcination 400–800 °C at 1–5 K/min','viscosity rises to gel point|capillary stress on drying sets cracking|dip/spin film 20–500 nm per coat',
'','drying shrinkage 50–90 % by volume|humidity 30–60 % RH',
'water:alkoxide ratio|pH/catalyst|ageing time|drying rate','particle size|porosity|surface area|crystalline phase',[6,7,8,9,10,21],'flammable alcohols|moisture-reactive alkoxides',
'Low-temperature route but solvent-rich; recover alcohol and avoid chlorinated precursors.',null,'gel time minutes–days');
e(4,'N 0 me ox ce lc L P I',[10,90],A,[60,7200],
['air or N₂ blanket','stirred aqueous solution','mechanical mixing; mild heating','bench to continuous multi-tonne plants'],
'salt 0.01–2 M|precipitant (NaOH, NH₃, carbonate) to supersaturation 10–10⁴|pH held ±0.2','near-isothermal, 10–90 °C|exotherm of neutralisation removed by jacket','micromixing time must be shorter than nucleation time|stirring 200–1500 rpm|ageing causes Ostwald ripening',
'','filtration, washing and drying downstream|agglomeration during drying',
'supersaturation|pH|mixing intensity|temperature','particle size|size distribution|phase|agglomeration',[6,7,8,9,21,26],'caustic reagents|fine powder handling',
'Simple and scalable; salt-laden wash water is the main burden — recycle mother liquor and recover by-product salts.',null,'kg–t per day (continuous)');
e(5,'N 0 ox me ce lc L P I',[20,90],A,[600,7200],
['N₂ or Ar blanket for Fe²⁺ systems','stirred aqueous solution','mechanical mixing; mild heating','bench to tonne-scale stirred tanks'],
'mixed cations at target stoichiometry (Fe³⁺:Fe²⁺ = 2:1 for Fe₃O₄)|pH 9–11 with NH₄OH or NaOH|total metal 0.05–1 M','isothermal 20–90 °C|optional calcination 300–900 °C for mixed oxides','rapid base addition gives small particles|vigorous stirring 500–1500 rpm|dissolved O₂ must be excluded for magnetite',
'magnetic decantation of product','washing to neutral pH|surface coating to prevent aggregation',
'cation ratio|pH|temperature|base addition rate','composition|particle size 5–20 nm|magnetisation|phase purity',[6,7,8,13,21,26],'ammonia vapour|caustic solutions',
'Aqueous and energy-light; treat ammonium/chloride effluent and prefer iron-based over cobalt- or nickel-rich compositions.',null,'g–kg per batch');
e(6,'N 0 1 me lc hp L P',[0,100],A,[60,14400],
['air (Au, Ag) or inert gas (Cu)','stirred aqueous or polyol solution','stirring; hotplate or ice bath','mL flasks to 100 L reactors'],
'metal salt 0.1–10 mM (HAuCl₄, AgNO₃)|reductant: citrate, NaBH₄, ascorbic acid, polyol|stabiliser: citrate, PVP, CTAB','boiling 100 °C (citrate) or 0–25 °C (borohydride)|polyol route 120–200 °C','burst nucleation needs fast mixing|seeded growth by slow addition',
'','glassware cleanliness critical (aqua regia)|storage dark, 4 °C',
'reductant strength and ratio|stabiliser|temperature|pH','diameter 2–100 nm|shape|plasmon wavelength|dispersity',[6,7,8,12,21],'borohydride releases H₂|silver and gold salts are corrosive/staining',
'Precious-metal use dominates cost and impact; recover metal from waste and favour benign reductants (citrate, ascorbate).',null,'mg–g per batch');
e(7,'N 2 me co cf lc L P I',[20,95],A,[60,7200],
['air','aqueous plating bath','chemical (autocatalytic reduction); bath heating','beakers to m³ plating tanks'],
'metal salt 5–30 g/L (NiSO₄, CuSO₄, Ag, Au)|reductant: hypophosphite, formaldehyde, DMAB|pH 4–5 (acid Ni–P) or 9–13|complexants and stabilisers','Ni–P 85–95 °C|Cu 25–60 °C|±1 K control for uniform rate','agitation and filtration|H₂ evolution at the surface|bath loading 0.5–2.5 dm²/L',
'none — no external current','surface activation (Pd/Sn) required on insulators|adhesion depends on pretreatment',
'bath temperature|pH|reductant concentration|stabiliser level','thickness|P or B content|uniformity|hardness',[6,7,8,11,21,26],'formaldehyde and nickel salts are toxic|hydrogen evolution',
'Uniform coatings without power, but baths have finite life; extend by replenishment and recover nickel from spent bath.',null,'2–25 µm/h');
e(8,'N 1 2 3 me ox sc al lc L P I',[15,80],A,[10,36000],
['air or N₂-purged electrolyte','liquid electrolyte, three-electrode cell','DC, pulsed or potentiostatic electrical power','cm² lab cells to m² plating lines'],
'metal ion 0.01–1 M|supporting electrolyte 0.1–1 M|pH 1–12 by system|additives ppm level','bath 15–80 °C, ±1 K|Joule heating small','diffusion layer 10–200 µm|rotating disc or agitation sets limiting current|template pore filling',
'current density 1–100 mA/cm²|potential −0.2 to −1.5 V vs reference|pulse 1–1000 Hz','substrate must be conductive|internal stress in deposits',
'potential/current density|bath composition|temperature|additives','thickness or wire length|grain size|composition|morphology',[6,7,9,11,16,21],'acid or cyanide electrolytes|hydrogen evolution',
'Efficient room-temperature route; manage metal-bearing rinse water with closed-loop recovery.',null,'0.01–1 µm/min');
e(9,'N 2 3 ce ox ca co cf lc L P I',[15,40],A,[10,1800],
['air','colloidal suspension in ethanol, isopropanol or water','DC or pulsed electric field','cm² coupons to automotive-scale coating lines'],
'solids 0.1–20 g/L|zeta potential |ζ| > 30 mV|charging additives (I₂, Mg salts, polyelectrolyte)|conductivity 1–100 µS/cm','ambient|post-sintering 400–1400 °C for ceramics','Hamaker law: mass ∝ mobility × field × concentration × time|water electrolysis causes bubbles above ~3 V',
'field 10–300 V/cm (organic) or < 5 V/cm DC (aqueous)|electrode gap 5–20 mm','conductive substrate|drying cracks above critical thickness',
'field strength|time|solids loading|zeta potential','coating thickness 0.1–100 µm|packing density|uniformity',[8,9,11,14,21,26],'flammable solvent with high voltage|nanoparticle suspensions',
'High material utilisation; prefer aqueous or ethanol suspensions and pulsed fields to limit electrolysis.',null,'0.1–10 µm/min');
e(10,'N 0 me ox sc hp L',[15,70],A,[600,86400],
['air or N₂','water-in-oil microemulsion (nanodroplet reactors)','stirring; mild heating','lab to small pilot batches'],
'water:surfactant molar ratio w₀ = 5–20 sets droplet size|surfactant AOT, CTAB, Triton X-100|oil: isooctane, cyclohexane|reactant 0.01–0.5 M in water pools','isothermal 20–50 °C within the microemulsion stability window','droplet collision and exchange ~µs–ms|two-microemulsion mixing',
'','product recovery by destabilisation and washing',
'w₀ ratio|reactant concentration|surfactant type|temperature','particle size 2–50 nm|narrow size distribution|shell thickness',[6,7,8,17,21],'flammable hydrocarbon oil|surfactant waste',
'Excellent size control but low yield per litre of solvent; recycle oil and surfactant, reserve for high-value particles.',null,'mg–g per batch');
e(11,'N 0 1 2 sc me ox hp L P',[20,360],A,[60,14400],
['Ar or N₂ on a Schlenk line','high-boiling organic solvent (octadecene, oleylamine) or water','heating mantle; rapid precursor injection','25 mL flasks to 100 L-class and flow reactors'],
'precursor 1–100 mM|ligands: oleic acid, oleylamine, TOPO, thiols|ligand:metal ratio 1–20','hot injection 240–320 °C then growth 200–300 °C|heat-up ramps 5–30 K/min|quench by solvent or bath','burst nucleation on injection, then diffusion-controlled size focusing|stirring 300–1200 rpm',
'in-situ absorbance/PL monitoring of size','air-free handling|purification by antisolvent precipitation',
'temperature|growth time|ligand chemistry|precursor reactivity','diameter 2–20 nm|size dispersion < 5–10 %|shape|emission wavelength',[6,7,8,12,17,21],'pyrophoric or toxic precursors (Cd, Pb, phosphines)|hot flammable solvent',
'Solvent- and ligand-intensive; move to heavy-metal-free compositions (InP, CuInS₂, ZnSe) and recycle solvent.',null,'g per batch; kg/day in flow');
e(12,'N 0 3 po co lc L P I',[20,150],[1e5,1e6],[1800,86400],
['N₂ blanket (radical routes are O₂-inhibited)','emulsion, miniemulsion, suspension or solution','heating/cooling jacket; thermal, redox or photo initiation','flasks to > 100 m³ industrial reactors'],
'monomer 5–50 wt%|initiator 0.1–2 wt% of monomer|surfactant above CMC|chain-transfer or RAFT/ATRP agent','emulsion polymerisation 50–90 °C|strongly exothermic (50–100 kJ/mol) — jacket cooling','viscosity rise and gel effect|stirring for heat removal|semi-batch monomer feed',
'','latex shear stability|coagulum control',
'surfactant concentration|initiator level|temperature|monomer feed rate','particle diameter 30–500 nm|molar mass|composition|filler dispersion',[6,7,8,9,10,26],'thermal runaway|toxic volatile monomers',
'Water-borne emulsion routes avoid solvent; strip residual monomer and consider bio-based monomers.',null,'conversion > 95 % in 1–8 h');
e(13,'N 2 po co cf lc L P I',[15,40],A,[5,600],
['air','interface between immiscible aqueous and organic phases','chemical (fast condensation reaction)','lab coupons to roll-to-roll membrane lines'],
'aqueous amine 0.1–3 wt% (m-phenylenediamine, piperazine)|organic acyl chloride 0.05–0.3 wt% (TMC in hexane)|acid acceptor','ambient; cure 60–90 °C for 2–10 min','diffusion of amine into the organic phase is self-limiting|film forms in seconds',
'','porous support (polysulfone) wetted with amine|thin film 10–200 nm is fragile unsupported',
'monomer concentrations|contact time|support pore structure|curing temperature','film thickness|cross-link density|permeance|salt rejection',[7,8,9,21,26],'acyl chlorides are corrosive|hexane is flammable',
'Low energy and fast; replace hexane with safer isoparaffins and capture solvent vapour.',null,'film forms in 5–60 s');
e(14,'N M 2 po co bi cf lc L P',[15,40],A,[600,86400],
['air','alternating aqueous polyelectrolyte/nanoparticle solutions and rinses','ambient; dip, spray or spin robots','coupons to spray-coated m² panels'],
'polyelectrolyte 0.1–5 g/L (PAH, PSS, PDDA)|ionic strength 0–1 M NaCl|pH sets charge density of weak polyelectrolytes','ambient','adsorption 1–20 min per layer by dipping, seconds by spraying|rinse between layers',
'','charged or primed substrate|humidity swells films',
'number of bilayers|ionic strength|pH|adsorption time','thickness (1–5 nm per bilayer)|composition profile|roughness|refractive index',[8,12,17,21],'low hazard; dilute polyelectrolytes',
'Water-based and ambient; rinse-water volume is the main burden — use spray assembly and recirculation.',null,'1–5 nm per bilayer');
e(15,'N 1 2 sc ox al cf lc L P I',[25,95],A,[600,43200],
['air, open or loosely covered bath','dilute aqueous solution','water bath or hotplate','beakers to m² solar-module baths'],
'metal salt 1–100 mM (Zn, Cd)|complexant NH₃, HMTA, citrate|chalcogen source thiourea for sulphides|pH 9–12 (or 5–7 with HMTA)','ZnO nanorods 60–95 °C for 1–6 h|CdS 60–80 °C for 5–30 min','heterogeneous growth competes with bulk precipitation|stirring improves uniformity|bath depletes — single use',
'','seeded or activated substrate|vertical substrate avoids sediment',
'temperature|pH/complexant ratio|precursor concentration|time','film thickness|rod length and diameter|density|orientation',[6,7,8,9,21,26],'cadmium and thiourea toxicity|ammonia vapour',
'Very low capital and energy; low precursor utilisation — reuse baths and replace Cd buffers with Zn(O,S).',null,'1–20 nm/min');
e(16,'N 0 1 me ox sc hp L',[20,250],A,[600,86400],
['air or inert gas; negligible solvent vapour pressure','ionic liquid (e.g. imidazolium salts)','conventional or microwave heating; electrodeposition','mL to L scale'],
'precursor 1–100 mM|ionic liquid acts as solvent, stabiliser and template|water content < 100 ppm for sensitive systems','stable to 250–400 °C|efficient microwave absorption','high viscosity 20–500 mPa·s slows diffusion|ionic conductivity 0.1–2 S/m',
'electrochemical window 4–6 V enables reactive-metal deposition','hygroscopic — dry handling|product washing is slow',
'anion/cation choice|temperature|water content|precursor concentration','particle size|shape|phase|dispersion stability',[6,7,8,11,17,21],'some ionic liquids are toxic or release HF (PF₆⁻, BF₄⁻)',
'Non-volatile but costly and not automatically green; recycle the ionic liquid and choose low-toxicity anions.',null,'mg–g per batch');
e(17,'N 0 1 me ox lc L P',[20,150],A,[600,86400],
['air','deep eutectic solvent (e.g. choline chloride:urea 1:2)','hotplate, microwave or electrodeposition','mL to tens of litres'],
'salt:hydrogen-bond-donor ratio 1:2|precursor 1–200 mM|water 0–20 wt% tunes viscosity','25–150 °C; urea-based DES decompose above ~100 °C','viscosity 50–1000 mPa·s at 25 °C|slow diffusion favours small particles',
'electrochemical window 2–3 V','hygroscopic|washing with water or ethanol',
'DES composition|water content|temperature|precursor concentration','particle size|shape|phase|surface chemistry',[6,7,8,11,21,28],'low acute hazard|ammonia from urea decomposition',
'Cheap, biodegradable components; viscosity raises mixing energy — recycle the solvent after product separation.',null,'mg–g per batch');
e(18,'N 0 3 ox ce po hp L P I',[31,450],[7.4e6,4e7],[1,14400],
['supercritical water, CO₂ or alcohol','supercritical fluid','preheaters and high-pressure pumps','mL batch cells to continuous t/yr plants'],
'metal salt 0.01–1 M|water critical point 374 °C, 22.1 MPa|CO₂ critical point 31.1 °C, 7.38 MPa','continuous hydrothermal 380–450 °C at 23–30 MPa|heating to reaction temperature in < 1 s by mixing','gas-like diffusivity, liquid-like density|residence time 1–60 s in flow|no surface tension',
'','pressure vessel codes apply|corrosion in supercritical water (Ni alloys)',
'temperature|pressure (density)|residence time|mixing geometry','particle size 5–100 nm|crystallinity|porosity|surface area',[6,7,9,10,21,26],'high pressure and temperature|corrosive supercritical water',
'No organic solvent and very fast kinetics; pumping and heating energy is high — use heat integration.',null,'residence 1–60 s; kg/h in flow');
e(19,'N 0 1 ox sc me lc L P',[60,250],[1e5,3e6],[60,3600],
['sealed vessel or open reflux','polar liquid (water, glycols, alcohols)','microwave dielectric heating','10 mL vials to L-scale batch and flow units'],
'precursor 1–200 mM|solvent loss tangent > 0.1 for efficient heating|same chemistry as hydro/solvothermal','volumetric heating 10–100 K/min|hold 1–60 min|IR or fibre-optic temperature sensing','penetration depth 1–3 cm in water limits vessel diameter|stirring evens hot spots',
'2.45 GHz, 100–1500 W|power density 1–50 W/mL','pressure-rated vessels with relief|metal parts cause arcing',
'temperature|hold time|power/ramp rate|solvent','particle size|aspect ratio|phase|reaction yield',[6,7,10,12,21],'rapid pressure rise|superheating',
'Cuts reaction time and energy by an order of magnitude at small scale; scale-up needs flow reactors.',null,'reactions complete in 1–60 min');
e(20,'N 0 me ox co lc L P',[0,80],A,[300,14400],
['air or Ar (Ar raises cavitation temperature)','liquid with acoustic cavitation','ultrasound horn or bath','mL to tens of litres; flow cells'],
'precursor 1–100 mM|volatile organometallics decompose inside bubbles|radicals H· and OH· from water sonolysis','bulk held 0–80 °C by cooling|collapsing bubbles ~5000 K and ~100 MPa for < 1 µs|cooling rates > 10⁹ K/s','microjets and shock waves break agglomerates|acoustic streaming mixes',
'20 kHz–1 MHz|intensity 10–100 W/cm² at the horn','horn tip erosion contaminates with Ti|noise enclosure',
'acoustic power|frequency|bulk temperature|dissolved gas','particle size|amorphous vs crystalline|dispersion|coating uniformity',[7,8,9,10,21],'noise|aerosol generation',
'Ambient bulk conditions and short times; acoustic efficiency is low, so use pulsed operation.',null,'mg–g per hour');
e(21,'N 0 2 sc me ox L',[10,60],A,[60,14400],
['air or inert, UV-transparent vessel','liquid solution or thin film','UV/visible lamp, LED or laser','cuvettes to falling-film and flow photoreactors'],
'precursor 0.1–10 mM|photoinitiator or sensitiser|radical scavengers and capping agents','near ambient; lamp heat removed by cooling jacket','Beer–Lambert attenuation limits optical path to mm–cm|thin-film or flow geometry for uniform dose',
'wavelength 254–450 nm|irradiance 1–100 mW/cm²|dose = irradiance × time','quartz or borosilicate windows|window fouling',
'wavelength|irradiance|exposure time|sensitiser concentration','particle size|shape (plasmon-mediated)|conversion|spatial pattern',[7,8,12,21],'UV exposure|ozone from short-wavelength lamps',
'Room-temperature and reagent-light; use efficient LEDs matched to the absorption band.',null,'mg–g per hour');
e(22,'N 0 me ox sc po hp L P I',[20,300],[1e5,1e7],[1,1800],
['closed channels, inert-purged feeds','single-phase or segmented liquid flow','heated chips or coils; syringe/HPLC pumps','µL/min chips to numbered-up kg/day systems'],
'same chemistries as batch at 1–500 mM|stoichiometry set by flow-rate ratio','heating/cooling within ms–s|isothermal to ±1 K|back-pressure regulator allows superheating','channels 0.1–1 mm, Re 1–500 (laminar)|mixing time 1–100 ms|segmented flow narrows residence-time distribution|fouling and clogging risk',
'inline UV–vis/PL feedback','pressure rating of tubing and fittings|pulsation-free pumping',
'residence time|temperature|flow-rate ratio|mixer type','particle size|dispersity|composition|reproducibility',[7,8,9,10,20,26],'small hold-up lowers hazard|leaks under pressure',
'Low inventory, high selectivity and easy heat recovery; scale by numbering-up rather than larger batches.',null,'residence 1 s–30 min');
e(23,'N 0 3 ox me ce co lc L P I',[20,150],A,[600,172800],
['sealed jar in air or Ar','dry or liquid-assisted solid state','mechanical impact and shear (ball mill, extruder)','g-scale jars to continuous twin-screw extrusion'],
'stoichiometric solid reactants|liquid-assisted grinding 0.1–1 µL/mg|salt diluent limits particle growth','jar bulk 20–150 °C|transient contact temperatures higher|self-propagating reactions possible','no bulk solvent transport; mixing by impact|ball-to-powder mass ratio 5:1–40:1',
'','milling media wear contaminates (Fe, ZrO₂, WC)|frequency 5–30 Hz or 200–600 rpm',
'milling energy/speed|time|ball-to-powder ratio|liquid additive','crystallite size 5–50 nm|phase|conversion|defect density',[6,7,14,16,21,26],'fine powder dust|exothermic self-propagating reactions',
'Solvent-free and near ambient; energy use per kg falls strongly with continuous extrusion.',null,'g/h (jars) to kg/h (extruders)');
e(24,'N 0 1 2 ce ox hp L P I',[300,1100],A,[1800,43200],
['air or inert gas over the melt','molten salt flux (chlorides, nitrates, sulphates, hydroxides)','resistance furnace','crucibles to industrial kilns'],
'salt:reactant mass ratio 1:1–20:1|NaCl–KCl eutectic 657 °C, LiCl–KCl 353 °C|reactant solubility in the melt controls mechanism','heating 2–10 K/min|hold 1–12 h|cooling rate affects crystal habit','dissolution–precipitation or template growth|fast ion diffusion (~10⁻⁹ m²/s) in the melt',
'','crucible corrosion (alumina, Pt)|salt removed by water washing',
'salt type|temperature|salt:reactant ratio|time','crystal size|shape (plates, rods, cubes)|phase purity|agglomeration',[6,7,10,16,21],'hot molten salt burns|chloride vapour corrosion',
'Lowers synthesis temperature by 200–400 K versus solid state; recover and reuse the salt from wash water.',null,'kg per batch');
e(25,'N 1 3 me ox sc ca po al pe lc L P',[15,100],A,[600,86400],
['air, liquid-filled template','pores of anodic alumina, track-etched membranes, mesoporous silica or colloidal crystals','electrodeposition, sol–gel infiltration, or vapour infiltration','cm² membranes to wafer scale'],
'template pore diameter 10–400 nm, length 1–100 µm|precursor infiltrated by electrodeposition, sol–gel, CVD/ALD|template removed by NaOH, HF or calcination','filling at 15–100 °C|calcination or template burn-out 400–600 °C','pore filling limited by diffusion and wetting|capillary infiltration|vacuum assists degassing of pores',
'electrodeposition 1–50 mA/cm² into pores','template brittleness|collapse of released wires on drying',
'pore diameter|pore length|filling time/charge|precursor concentration','wire diameter|length|array pitch|wall thickness of tubes',[4,7,8,9,11,21],'HF or caustic template etchants',
'Templates are consumed; anodic alumina is inexpensive, but minimise etchant use and prefer reusable imprinted moulds.',null,'fill rate 0.01–1 µm/min');

// ── B. Thermal and high-temperature environments ──
e(26,'N 3 0 ce ox co lc L P I',[300,1800],A,[1800,172800],
['air, or flowing N₂/Ar','solid-state powder or compact','resistance elements (Kanthal, SiC, MoSi₂)','lab muffle furnaces to tunnel kilns'],
'precursor decomposition (carbonates, nitrates, hydroxides)|stoichiometry fixed at mixing|oxygen partial pressure of furnace gas','ramp 1–20 K/min|hold 1–24 h at 0.5–0.8 of melting temperature|furnace cooling','solid-state diffusion controls reaction|gas release needs open porosity|bed depth affects uniformity',
'','crucible/setter compatibility|grain growth competes with densification',
'peak temperature|hold time|heating rate|atmosphere','crystallite size|phase|density|surface area',[6,10,14,16,21,26],'hot surfaces|NOₓ/CO₂ off-gas',
'Energy-intensive; use lower-temperature precursors, shorter cycles and electrified, well-insulated kilns.',null,'kg–t per batch');
e(27,'N M 2 1 ca sc ox me cf hp ep L P I',[300,1200],[10,1.013e5],[60,14400],
['precursor in H₂/Ar/N₂ carrier','gas phase reacting on a heated substrate or catalyst','hot-wall tube furnace or cold-wall susceptor','cm² samples to 300 mm wafers and roll-to-roll foil'],
'precursors: CH₄, C₂H₄, SiH₄, metal-organics, halides|graphene: CH₄:H₂ ≈ 1:10–1:1000 on Cu|CNT: Fe/Co/Ni catalyst 1–10 nm','graphene on Cu 1000–1050 °C|CNT 600–900 °C|ramp 10–50 K/min, fast cool','boundary-layer diffusion vs surface-reaction control|total flow 10–2000 sccm|precursor depletion along the tube',
'','substrate thermal expansion mismatch on cooling|catalyst foil surface quality',
'temperature|precursor partial pressure|pressure|time','thickness/layer number|grain size|tube diameter and length|defect density',[6,7,9,10,16,17,21],'flammable/pyrophoric gases (H₂, SiH₄)|toxic by-products',
'High-purity films at high thermal cost; reduce temperature with plasma or catalysts and abate exhaust gases.',null,'0.1–100 nm/min');
e(28,'N M 2 1 sc ca ce cf L P I',[100,700],[10,1300],[60,7200],
['reactive gases with Ar/He dilution','low-pressure glow-discharge plasma','RF or microwave plasma plus substrate heater','wafer tools to m² display and solar lines'],
'SiH₄ + NH₃/N₂O for SiNₓ/SiO₂|CH₄/C₂H₂ + H₂ for carbon|radical and ion chemistry replaces thermal activation','substrate 200–400 °C (dielectrics), 500–700 °C (CNT, graphene)|plasma heating of the surface','showerhead flow 100–5000 sccm|residence time 0.1–1 s|powder formation at high power',
'RF 13.56 MHz (and 40–400 kHz)|power density 0.01–1 W/cm²|sheath bias aligns CNTs vertically','film stress −1 to +1 GPa tuned by frequency mixing|hydrogen content affects stability',
'substrate temperature|RF power|gas ratio|pressure','thickness|stress|refractive index|hydrogen content',[7,10,11,12,16,21],'silane is pyrophoric|toxic and greenhouse process gases',
'Enables low-temperature deposition; replace NF₃/SF₆ chamber cleans or abate them (high global-warming potential).',null,'10–500 nm/min');
e(29,'N M 2 ox ce me sc cf hp L P I',[50,400],[10,1000],[600,36000],
['alternating precursor pulses and N₂/Ar purges','gas–surface self-limiting half-reactions','heated chamber; optional plasma','single wafers to batch and spatial roll-to-roll ALD'],
'precursors: TMA, DEZ, TiCl₄, TDMAT, metal amidinates|co-reactant H₂O, O₃, O₂ plasma, NH₃|doses exceed saturation','ALD window 50–400 °C; Al₂O₃ typically 150–300 °C|isothermal to ±2 K','pulse 0.01–1 s, purge 1–30 s|conformal in aspect ratios > 100:1 given longer exposure|Knudsen diffusion in pores',
'optional remote plasma 100–300 W','nucleation delay on inert surfaces|pinhole-free barrier films',
'number of cycles|temperature|precursor/co-reactant|pulse and purge times','thickness (0.05–0.15 nm per cycle)|conformality|composition|crystallinity',[7,8,9,17,21,26],'pyrophoric metal alkyls|corrosive by-products (HCl)',
'Low precursor utilisation (often < 30 %) and slow; spatial ALD and precursor recovery improve both.',null,'0.05–0.15 nm/cycle; 0.3–6 nm/min');
e(30,'N M 2 me ox ce sc hp pe L P I',[20,600],[1e-5,5],[60,14400],
['high vacuum, with Ar or reactive gas when sputtering','vapour from a solid source condensing on a substrate','resistive, electron-beam, arc or sputter source','R&D coaters to architectural-glass lines'],
'source purity 99.9–99.999 %|reactive gas (O₂, N₂) partial pressure sets stoichiometry|base pressure < 10⁻⁴ Pa limits contamination','substrate 20–600 °C; structure zones follow T/T_melt|radiant heating from the source','line-of-sight transport; mean free path > chamber size below ~10⁻² Pa|planetary rotation for uniformity',
'substrate bias 0 to −200 V|ion assistance densifies films','film stress and adhesion|shadowing gives columnar growth',
'deposition rate|substrate temperature|working pressure|bias','thickness|density|stress|grain structure',[7,10,12,16,17,21],'high voltage|fine metal dust during maintenance',
'Solvent-free and material-efficient with shields recycled; pump-down and heating energy dominate.',null,'0.1–10 nm/s');
e(31,'N M 2 me sc po hp lc L P I',[20,300],[1e-5,1e-3],[60,3600],
['high vacuum','vapour from a resistively heated boat or crucible','Joule-heated W/Mo/Ta boat or Knudsen cell','bell jars to roll-to-roll web metallisers'],
'source metals Al, Ag, Au, Cr, organics|alloys fractionate by vapour pressure','source 800–2000 °C to reach ~1 Pa vapour pressure|substrate near ambient to 300 °C','line-of-sight, cosine-law flux|throw distance 20–50 cm|low adatom energy ~0.1 eV',
'boat current 50–400 A at a few volts','poor step coverage — suits lift-off|adhesion layer (Cr, Ti) 2–5 nm',
'source power/rate|thickness (quartz monitor)|substrate temperature|geometry','thickness|grain size|roughness|purity',[7,10,12,17,21],'hot sources|metal flakes',
'Simple and solvent-free; low material utilisation on small substrates — recover precious metals from shields.',null,'0.1–5 nm/s');
e(32,'N M 2 me ox ce sc hp pe L P I',[20,400],[1e-6,1e-3],[60,7200],
['high vacuum; optional O₂ backfill ~10⁻² Pa','vapour from an electron-beam-heated melt','electron gun 5–10 kV, 1–15 kW','R&D chambers to large optical-coating plants'],
'refractory metals and oxides (Ti, Pt, SiO₂, TiO₂, HfO₂)|oxides lose oxygen — reactive backfill','local melt > 2000 °C in a water-cooled hearth|substrate 20–400 °C','line-of-sight flux|planetary fixtures; thickness uniformity ±1–2 %|ion-assisted deposition densifies',
'beam sweep pattern|ion source 50–300 eV for assistance|X-ray emission from the gun','porous columnar films without ion assist|humidity shifts optical thickness',
'beam power/rate|substrate temperature|O₂ backfill|ion-assist energy','thickness (±1 nm)|refractive index|packing density|stress',[7,10,12,16,21],'high voltage and X-rays|hot hearth',
'High-purity, high-rate route; energy per wafer is modest once under vacuum — maximise batch loading.',null,'0.05–10 nm/s');
e(33,'N M 2 me ox ce sc co pe hp L P I',[20,600],[0.1,5],[60,14400],
['Ar, with O₂/N₂ for reactive sputtering; base < 10⁻⁴ Pa','magnetically confined glow-discharge plasma','DC, pulsed-DC, RF or HiPIMS power','R&D tools to in-line coaters for glass and web'],
'target purity 99.9–99.999 %|reactive gas flow shows hysteresis (poisoned vs metal mode)|alloy composition transfers closely','substrate 20–600 °C|plasma heat load 0.1–1 W/cm²|target water-cooled','mean free path ~mm–cm at 0.1–5 Pa|target–substrate distance 5–15 cm|sputtered-atom energy 1–10 eV',
'DC 300–700 V|power density 1–20 W/cm²|RF 13.56 MHz for insulators|magnetic field 20–50 mT at the target','compressive or tensile stress tuned by pressure and bias|good adhesion',
'power|working pressure|reactive gas flow|substrate bias','thickness|stoichiometry|stress|density and grain size',[7,10,11,12,16,21],'high voltage|toxic targets (Be, Cd, Pb compounds)',
'Scalable and solvent-free; target utilisation is 30–40 % for planar and ~80 % for rotatable targets — recycle spent targets.',null,'0.1–2 nm/s');
e(34,'N M 2 0 sc ox me ep hp L P',[200,800],[1e-9,1e-4],[1800,86400],
['ultra-high vacuum; base ~10⁻⁸ Pa','molecular beams on a heated single-crystal substrate','effusion cells, crackers, e-beam sources; RHEED monitoring','single-wafer R&D to multi-wafer production'],
'elemental sources 6N–7N purity|V/III beam-equivalent pressure ratio 5–50|dopant cells (Si, Be)','GaAs 550–620 °C, Si/SiGe 400–700 °C, oxides 600–900 °C|substrate ±1 K','collision-free molecular flow|shutter switching < 0.5 s gives monolayer abruptness|surface migration length sets morphology',
'RHEED 10–30 keV for in-situ growth-rate oscillations','lattice mismatch strain drives quantum-dot formation (Stranski–Krastanov)|LN₂ cryoshrouds',
'substrate temperature|flux ratios|growth rate|shutter sequence','layer thickness (monolayer)|composition|doping|interface abruptness',[3,7,11,12,16,17,21],'toxic As/P sources|pyrophoric phosphorus residues',
'Unmatched purity at high energy and LN₂ cost; long campaigns and closed-cycle cryocooling reduce the footprint.',null,'0.1–1 µm/h (≈ 0.1–1 monolayer/s)');
e(35,'N M 2 ox ce co sc ep hp L',[20,900],[1e-6,50],[300,14400],
['vacuum or O₂/N₂/Ar background 0.1–50 Pa','laser-ablation plasma plume onto a heated substrate','pulsed UV excimer or Nd:YAG laser','cm² research samples; 200 mm tools exist'],
'stoichiometric ceramic target|oxygen background controls oxidation|multi-target sequencing for superlattices','substrate 500–850 °C for epitaxial oxides|post-anneal in 10⁴–10⁵ Pa O₂','plume expansion forward-peaked (cosⁿ, n ≈ 5–20)|target–substrate 4–8 cm|particulates (droplets) a known defect',
'KrF 248 nm, 10–30 ns|fluence 1–5 J/cm²|repetition 1–20 Hz','thermal-mismatch and epitaxial strain|small uniform area unless scanned',
'fluence|background pressure|substrate temperature|pulse count','thickness|stoichiometry|epitaxial quality|oxygen content',[3,7,10,12,16,17],'class 4 UV laser|toxic excimer gases',
'Excellent stoichiometry transfer for complex oxides at low throughput; best suited to research and niche production.',null,'0.01–0.1 nm/pulse');
e(36,'N 0 ox ca ce lc P I L',[1500,2500],A,[0.001,1],
['combustion gases (fuel + O₂/air)','flame aerosol','combustion enthalpy','lab burners (g/h) to t/h industrial reactors'],
'precursors SiCl₄, TiCl₄, HMDSO, metal-organic sprays|equivalence ratio sets oxidising/reducing conditions|precursor loading 0.1–1 M (spray)','peak 1500–2500 °C|cooling 10⁴–10⁶ K/s|temperature–time history sets sintering','residence time 1–1000 ms|coagulation and sintering in free-molecular regime|turbulent mixing',
'optional external electric field narrows size','filter/baghouse collection|agglomerate strength',
'precursor feed rate|fuel/oxidant ratio|flame temperature|residence time','primary particle size 5–100 nm|aggregate size|crystal phase (anatase/rutile)|surface area',[6,7,9,10,21,26],'nanoparticle aerosol exposure|HCl from chloride precursors',
'Highest-volume nanoparticle route (carbon black, fumed silica, TiO₂); recover heat and scrub HCl/Cl₂ for recycle.',null,'g/h (lab) to t/h');
e(37,'N 0 ox ce co me P I L',[400,1500],A,[0.1,60],
['carrier air, N₂ or forming gas','aerosol droplets in a hot-wall reactor','tubular furnace or flame','lab tubes to industrial towers'],
'salt solution 0.01–2 M (nitrates, acetates, chlorides)|one droplet → one particle, composition preserved','evaporation, precipitation, decomposition, sintering in sequence|heating 10–1000 K/s','droplets 1–10 µm from ultrasonic or two-fluid nebulisers|residence 1–10 s|hollow particles when solute precipitates at the surface',
'ultrasonic nebuliser 1.7–2.4 MHz','collection by filter, cyclone or electrostatic precipitator',
'solution concentration|droplet size|temperature profile|residence time','particle size 0.05–2 µm|hollow vs dense|composition|crystallinity',[6,7,9,10,21,26],'NOₓ off-gas|aerosol exposure',
'Continuous and compositionally exact; evaporating water costs energy — use concentrated feeds and heat recovery.',null,'g/h to 100 kg/h');
e(38,'N 2 0 ox sc ce lc cf L P I',[250,1000],A,[1,3600],
['air or N₂ carrier','spray droplets onto a hot substrate or into a furnace','hotplate or furnace; pneumatic or ultrasonic atomiser','lab hotplates to in-line glass coaters'],
'precursor 0.01–0.5 M (chlorides, acetates, acetylacetonates)|solvent water/alcohol|dopants added in solution','substrate 250–600 °C for films|droplet must evaporate just before impact (CVD-like regime)','spray rate 1–20 mL/min|nozzle distance 15–40 cm|substrate cooling by spray',
'','thermal shock of glass|film adhesion and haze',
'substrate temperature|precursor concentration|spray rate|nozzle distance','film thickness|grain size|resistivity/transparency|roughness',[7,9,10,12,21,26],'HCl or solvent vapour|aerosol',
'Inexpensive open-air coating (e.g. FTO glass); capture acid off-gas and overspray.',null,'10–200 nm/min');
e(39,'N 0 sc ce ca hp L P',[800,2000],[1e4,1e5],[0.001,0.1],
['Ar/He with gaseous precursors','laser-heated gas jet, no hot walls','CW CO₂ laser 10.6 µm, 0.1–5 kW','g/h lab to kg/h pilot'],
'precursors SiH₄, C₂H₄, NH₃, Fe(CO)₅, TTIP|sensitiser (SF₆, C₂H₄) when precursor does not absorb at 10.6 µm','reaction zone 800–2000 °C|heating and cooling ~10⁵–10⁶ K/s','residence 1–50 ms in the beam|confinement by coaxial inert sheath flow',
'laser intensity 10²–10⁴ W/cm²','wall-free reaction gives high purity|filter collection',
'laser power|precursor flow|pressure|sensitiser fraction','particle size 5–50 nm|composition|crystallinity|size distribution',[7,9,10,12,21],'pyrophoric silane and nanopowders|toxic carbonyls',
'Clean, narrow-size powders at modest efficiency (laser wall-plug ~10 %); best for high-value Si, SiC, Si₃N₄.',null,'10 g/h–kg/h');
e(40,'N 1 0 ca me L P',[3000,6000],[1e4,1e5],[60,3600],
['He or Ar (H₂, N₂ additions); also in liquid','arc plasma between consumable electrodes','DC arc 50–150 A at 20–35 V','lab chambers; g/h'],
'graphite anode with Fe/Co/Ni/Y catalyst 1–5 at% for single-wall tubes|He 13–66 kPa','plasma 3000–6000 °C|steep gradients to water-cooled walls give rapid quench','anode sublimes, vapour condenses on cathode and walls|convection sets quench rate',
'current 50–150 A|electrode gap 1–3 mm held by feed','product is soot needing purification (oxidation, acid)|electrode consumption',
'arc current|gas type and pressure|catalyst composition|electrode gap','tube diameter|wall number|yield|crystallinity',[6,7,10,11,21],'intense UV and heat|carbon nanotube and metal fume inhalation',
'Highly crystalline tubes but energy-hungry and low-yield; purification acids add waste.',null,'0.1–10 g/h');
e(41,'N 0 ce me ox ca hp P I L',[3000,12000],[2e4,1.1e5],[0.001,0.1],
['Ar, N₂, H₂, He or reactive mixtures','thermal (equilibrium) plasma jet','DC arc torch or RF induction torch, 10 kW–MW','pilot to industrial kg/h'],
'solid, liquid or gas feed fully vaporised|reactive gas (N₂, CH₄, O₂) forms nitrides, carbides, oxides','plasma core ~10 000 K|quench 10⁵–10⁶ K/s by gas injection|steep radial gradients','residence 1–100 ms|homogeneous nucleation in the quench zone|incomplete vaporisation of coarse feed',
'RF 0.2–4 MHz induction (electrodeless, high purity)','water-cooled reactor walls|passivation of pyrophoric powders',
'plasma power|feed rate|quench gas flow|reactive gas','particle size 10–200 nm|phase|sphericity|purity',[6,7,9,10,21,26],'intense UV, noise, high power|pyrophoric nanopowders',
'Very high specific energy (often > 10 kWh/kg); justified for refractory or high-purity powders.',null,'0.1–100 kg/h');
e(42,'N 0 me hp L',[20,20000],A,[60,36000],
['flowing N₂, Ar or He at ambient pressure','gas-phase aerosol from electrode erosion','repetitive capacitor spark discharge','mg/h lab generators; parallelised to g/h'],
'electrode material is the only precursor — very high purity|alloy/mixed particles from dissimilar electrodes','carrier gas near 20 °C; each microsecond spark reaches ~20 000 K (window spans both)|quench > 10⁷ K/s','carrier flow 1–50 L/min dilutes and limits coagulation|primary particles 1–10 nm',
'breakdown 1–10 kV|spark energy 1–200 mJ|frequency 1 Hz–1 kHz (to 25 kHz)','electrode gap 0.5–3 mm|deposition by impaction, diffusion or electrostatic precipitation',
'spark energy|spark frequency|carrier gas flow|gas type','primary particle size|agglomerate size|composition|mass rate',[7,9,10,11,21],'ozone/NOₓ if air is present|nanoparticle aerosol',
'No solvents, precursors or waste streams; throughput is small, so scale by multiplexing electrode pairs.',null,'0.1–100 mg/h per electrode pair');
e(43,'N 2 sc me ox hp L P I',[400,1200],[100,1.1e5],[1,300],
['N₂, Ar, O₂, forming gas or NH₃','solid film on a wafer, radiatively heated','tungsten–halogen or arc lamps; lasers for ms annealing','single-wafer tools to 300 mm production'],
'dopant activation, silicidation, oxidation, crystallisation|ambient O₂ controlled to ppm','ramp 20–250 K/s|soak 1–120 s or spike|cooling 20–100 K/s','cold-wall chamber; gas flow 1–20 slm|diffusion length √(Dt) kept to nm',
'lamp emission 0.4–4 µm; emissivity-corrected pyrometry','thermal stress and wafer slip above ~1000 °C|pattern-dependent heating',
'peak temperature|soak time|ramp rate|ambient gas','grain size|dopant activation|phase (e.g. silicide)|junction depth',[10,11,14,16,21],'hot wafers|toxic/flammable ambients',
'Thermal budget is delivered only to the wafer for seconds, saving energy versus furnace anneals.',null,'1–5 min per wafer');
e(44,'N M 3 ce me co ox P I L',[400,2000],[1,1.1e5],[300,43200],
['vacuum, Ar, N₂ or air (conventional)','powder compact; graphite die for SPS','pulsed DC through die and powder (SPS) or furnace heating','cm discs to 300+ mm SPS and industrial furnaces'],
'nanopowder purity and surface oxides|sintering aids 0.1–5 wt%|carbon pick-up from graphite tooling','SPS 50–600 K/min to ≤ ~2000 °C, hold 1–20 min|conventional 1–10 K/min, hold 1–10 h|radial gradients in large dies','pore closure by grain-boundary and lattice diffusion|outgassing during ramp',
'pulsed DC 1–10 kA at < 10 V','uniaxial pressure 30–100 MPa (graphite dies); to 1 GPa with special tooling|density > 98 % with grains < 100 nm achievable',
'temperature|pressure|heating rate|hold time','relative density|grain size|hardness|residual porosity',[10,11,14,16,21,22],'high current|hot pressurised tooling',
'Short cycles cut energy several-fold versus hot pressing; graphite tooling is consumable.',null,'cycle 10–60 min');
e(45,'M N 3 me ce co sc al pe P I L',[200,1700],[1e-2,1.1e5],[3600,172800],
['vacuum or Ar','melt solidifying along an imposed temperature gradient','resistance or induction furnace with chill plate (Bridgman)','lab rods to turbine-blade casting furnaces'],
'alloy or eutectic composition sets phase fractions|impurity segregation coefficient k','gradient G = 1–20 K/mm|withdrawal V = 1–500 µm/s|melt superheat 50–200 K above liquidus','solute boundary layer D/V|constitutional undercooling when G/V is low|buoyant convection causes freckles',
'','mould/crucible compatibility|eutectic spacing λ²V = constant (Jackson–Hunt)',
'withdrawal rate|temperature gradient|composition|convection','lamellar/rod spacing 0.1–10 µm|alignment|phase fraction|texture',[6,10,14,16,18,21],'molten metal|vacuum furnace hazards',
'Long high-temperature cycles; self-organised eutectics give aligned composites without lithography.',null,'1–500 µm/s');

// ── C. Electrical, electromagnetic and plasma environments ──
e(46,'N 1 2 ox sc me al lc L P',[60,250],[1e5,4e6],[1800,86400],
['sealed aqueous or solvent solution; vapour-phase variants exist','liquid growth solution between electrodes','oven heating plus DC, AC or pulsed electric field','lab cells to instrumented pilot reactors'],
'same chemistry as hydrothermal/bath growth (1–100 mM)|field drives ion migration (Zn²⁺, OH⁻) to the growth front|pH 6–12','isothermal 60–250 °C|Joule heating σE² must stay small','electromigration adds to diffusion (Nernst–Planck)|electro-osmotic and electrothermal flow|porous separators shape the flux',
'field 0.1–100 kV/m|DC, pulsed or AC 1 Hz–1 MHz|electrode gap 1–50 mm','insulated or inert electrodes avoid electrolysis above ~1.23 V|vessel feedthroughs rated for pressure',
'field strength and direction|temperature|precursor concentration|time','orientation/alignment|length|growth rate|areal density',[7,9,10,11,20,21],'pressurised hot vessel|electrical hazard and H₂/O₂ from electrolysis',
'Field assistance raises rate and alignment for a small electrical input; avoid electrolysis, which wastes charge and makes gas.',null,'growth 2–30 nm/min');
e(47,'N M 2 1 3 me co cf lc L P I',[15,70],A,[1,7200],
['air; filtered, agitated bath','aqueous or ionic-liquid electrolyte','DC, pulse or pulse-reverse rectifier','lab cells to wafer plating and reel-to-reel lines'],
'metal ion 0.05–1.5 M (Cu, Ni, Au, Ag, Co)|suppressor, accelerator, leveller additives|pH and chloride controlled','bath 20–60 °C ±1 K','limiting current from diffusion layer|superconformal filling of vias|hydrogen bubbles create pits',
'current density 1–500 mA/cm²|pulse on-time 0.1–10 ms|overpotential sets nucleation density','seed layer required|deposit stress and hydrogen embrittlement',
'current density/waveform|additives|bath composition|temperature','grain size (to < 50 nm)|thickness|twin density|composition',[6,7,9,11,16,21],'acid and cyanide baths|metal-bearing effluent',
'Additive and efficient (current efficiency often > 90 %); closed-loop rinse and metal recovery are standard practice.',null,'0.01–2 µm/min');
e(48,'N M 1 0 ca me sc al pe L',[15,40],A,[1,600],
['air','dilute suspension over microelectrodes','AC electric field','chip-scale; wafer-scale parallel electrode arrays'],
'low-conductivity medium 10⁻⁴–10⁻² S/m|surfactant-stabilised nanotubes or particles at ng–µg/mL','ambient; electrothermal heating of a few K','force ∝ r³ Re[K(ω)] ∇|E|²|AC electro-osmosis competes below ~10 kHz|capture radius a few µm',
'1–20 V peak-to-peak|1 kHz–10 MHz|field 10⁵–10⁷ V/m|Clausius–Mossotti factor sign sets attraction/repulsion','electrode gap 1–50 µm|capillary forces on drying can misalign',
'frequency|voltage|medium conductivity|time','alignment|number of bridged objects|position|metallic/semiconducting selectivity',[8,9,11,12,20,21],'low hazard; nanotube suspensions',
'Places pre-made nanostructures with negligible energy and material use.',null,'assembly in 1–300 s');
e(49,'N 0 2 po ox me bi lc L P',[15,400],A,[60,14400],
['air, N₂ or CO₂ sheath','charged droplets in a cone-jet spray','high-voltage DC supply; syringe pump','single needles to multi-nozzle arrays'],
'solution conductivity 10⁻⁴–10⁻¹ S/m|solute 0.01–5 wt%|solvent volatility sets particle morphology','ambient for particles; substrate 200–400 °C for electrostatic spray deposition','flow 0.05–5 mL/h|droplets 0.1–10 µm (scaling d ∝ Q^(1/3))|Coulomb fission as solvent evaporates',
'2–20 kV|needle–collector 2–20 cm|cone-jet mode window in voltage and flow','humidity affects evaporation|grounded collector',
'flow rate|voltage|conductivity|concentration','particle size|morphology|coating thickness|porosity',[8,9,11,21],'high voltage|solvent vapour',
'Near-total material utilisation at ambient conditions; throughput needs multiplexed emitters.',null,'0.05–5 mL/h per nozzle');
e(50,'N 1 3 po ox ce co ca al lc L P I',[15,40],A,[600,86400],
['air at 20–60 % relative humidity','polymer solution or melt jet','high-voltage DC; syringe pump','single needle to needleless industrial lines'],
'polymer 5–25 wt% above entanglement concentration|solvent DMF, ethanol, water, HFIP|conductivity raised by salts|sol–gel precursor added for ceramic fibres','ambient; calcination 400–800 °C converts to oxide fibres','flow 0.1–5 mL/h|whipping instability thins the jet|solvent evaporation in flight',
'5–30 kV|tip–collector 5–25 cm (≈ 1 kV/cm)|rotating drum or gap collectors align fibres','mats 10–500 µm thick, porosity 70–90 %|humidity controls surface pores',
'solution concentration/viscosity|voltage|flow rate|distance','fibre diameter 50–1000 nm|alignment|bead formation|mat porosity',[8,9,11,14,21,26],'high voltage|toxic solvent vapour (DMF, HFIP)',
'Simple and scalable; solvent is 75–95 % of the feed — use benign solvents with vapour recovery or melt electrospinning.',null,'0.1–1 g/h per needle; kg/h needleless');
e(51,'N M 2 po ca ox me cf L P I',[20,150],[10,1.1e5],[5,1800],
['O₂, Ar, N₂, NH₃, H₂ or fluorocarbons; low-pressure or atmospheric','non-thermal plasma contacting a surface','RF, microwave, DBD or corona discharge','benchtop cleaners to in-line web treaters'],
'radicals and ions graft –OH, –COOH, –NH₂ or –F groups|treatment depth a few nm|hydrophobic recovery over hours–days','gas 20–150 °C; electrons 1–10 eV','radical flux and exposure time set the dose|shadowed and porous regions treated less',
'RF 13.56 MHz, 10–300 W|DBD 1–20 kV at kHz|ion energy 10–100 eV','surface roughening by etching|substrate must tolerate vacuum or ozone',
'gas chemistry|power|time|pressure','surface energy/contact angle|functional group density|roughness|adhesion',[8,11,17,21,23],'ozone and UV|fluorinated gases',
'Dry replacement for wet chemical primers; milliwatt-hours per cm².',null,'treatment 5 s–30 min');
e(52,'N M 2 sc ce ca po cf L P I',[20,400],[0.5,300],[60,7200],
['Ar with reactive gases','capacitively coupled RF glow discharge','13.56 MHz RF generator with matching network','R&D chambers to large-area deposition and etch tools'],
'precursor fragmentation by electron impact|deposition vs etching set by gas (SiH₄, CH₄, CF₄, O₂)|plasma polymerisation of organic monomers','substrate 20–400 °C|gas near ambient, electrons 2–5 eV','flow 10–1000 sccm|residence 0.1–10 s|dust formation at high power and pressure',
'13.56 MHz|0.05–2 W/cm²|self-bias 50–500 V|plasma density 10⁹–10¹⁰ cm⁻³','ion bombardment controls film density and stress',
'RF power|pressure|gas mixture|substrate temperature','thickness|composition|stress|etch or deposition rate',[7,11,12,16,21],'RF exposure|toxic process gases',
'Low-temperature processing; abate fluorinated and silane exhaust.',null,'5–200 nm/min');
e(53,'N 2 3 0 ca hp ep L P I',[600,1200],[2e3,3e4],[3600,360000],
['H₂ with 0.5–5 % CH₄ (N₂, O₂, B additions)','microwave plasma ball above the substrate','2.45 GHz or 915 MHz magnetron, 1–100 kW','25 mm R&D to 150+ mm production reactors'],
'atomic hydrogen etches sp² carbon and stabilises diamond|CH₄/H₂ ratio sets grain size (nano- vs microcrystalline)|diamond seeding 10¹⁰–10¹² cm⁻²','substrate 700–1100 °C|gas temperature ~3000 K in the plasma|power density 5–200 W/cm³','diffusion of H and CH₃ across the boundary layer|flow 100–1000 sccm',
'2.45 GHz, 1–6 kW (lab)|resonant cavity mode (TM₀₁ family)','thermal stress on cooling|substrate holder cooling',
'methane fraction|pressure|microwave power|substrate temperature','grain size|growth rate|sp³/sp² ratio|doping and colour centres',[6,7,10,12,17,21],'hydrogen flammability|microwave leakage',
'Electricity-dominated; lower-carbon than HPHT when run on low-carbon power, but long runs demand high uptime.',null,'0.1–10 µm/h (to ~100 µm/h single crystal)');
e(54,'M N 2 3 sc ox me ce pe L P I',[-110,400],[0.1,10],[10,3600],
['Cl₂, BCl₃, SF₆, CF₄, CHF₃, HBr, O₂, Ar','high-density inductively coupled plasma','ICP coil plus independent substrate RF bias','single-wafer tools up to 300 mm'],
'halogen radicals form volatile etch products|sidewall passivation (polymer, oxide) gives anisotropy|selectivity to mask 2–100:1','chuck −110 °C (cryogenic Si) to 400 °C (InP)|He backside cooling','ion flux ~10¹⁶ cm⁻² s⁻¹|aspect-ratio-dependent etch lag|loading and microloading',
'ICP 100–3000 W at 2–13.56 MHz|bias 0–500 W; ion energy 20–500 eV|density 10¹¹–10¹² cm⁻³','electrostatic clamping|plasma-induced damage depth 1–20 nm',
'ICP power|bias power|gas chemistry|pressure','etch depth|sidewall angle|selectivity|surface roughness',[5,7,11,12,21,23],'toxic corrosive gases (Cl₂, BCl₃, HBr)|RF power',
'Dry and precise; fluorinated gases have very high global-warming potential — abate above 95 %.',null,'etch 0.05–5 µm/min; resolution set by the mask');
e(55,'N M 1 3 me ox co al pe lc L P',[15,80],A,[1,3600],
['air or liquid carrier','magnetic particles in liquid, gel or curing resin','permanent magnets or electromagnets','mm samples to roll-to-roll aligned composites'],
'superparamagnetic Fe₃O₄ or ferromagnetic Ni/Co particles|steric or electrostatic stabilisation|matrix cure locks the structure','ambient; thermal energy kT competes with dipolar energy','chain growth by diffusion-limited aggregation|viscous drag sets assembly time|sedimentation',
'flux density 1–500 mT|gradient 1–100 T/m for transport|dipolar coupling parameter λ > 1 needed for chains','matrix viscosity window before gelation|chains give anisotropic stiffness and conductivity',
'field strength|field direction/gradient|particle concentration|viscosity','chain length|spacing|alignment|anisotropy',[8,9,13,14,20,21],'strong magnets (pinch, implants)',
'Contactless, low-energy alignment; iron-oxide particles are abundant and benign.',null,'chains form in 0.1–100 s');
e(56,'M N 3 2 co po me ox al pe L',[15,60],A,[1,3600],
['air or sealed liquid cell','colloidal suspension under combined electric, magnetic or optical fields','AC electrodes, coils or lasers','chip-scale cells; cm² with patterned electrodes'],
'particle volume fraction 0.1–30 %|medium permittivity and conductivity contrast|index-matched solvents for optical access','ambient; absorption heating limits optical power','dipolar chaining then coarsening into body-centred-tetragonal lattices|field-induced flows',
'electric field 10³–10⁶ V/m at 100 Hz–10 MHz|magnetic 1–100 mT, static or rotating|optical trapping 1–100 mW per trap','structures relax when the field is removed unless fixed by gelation or UV cure',
'field amplitude|frequency|field geometry|volume fraction','lattice symmetry|spacing|reconfiguration time|optical/EM response',[5,8,12,13,20,21],'laser and high-voltage hazards',
'Reversible, reprogrammable assembly at negligible material cost; areas are small.',null,'reconfiguration ms–minutes');
e(57,'M N 2 3 me po co lc L P',[15,80],A,[60,14400],
['air','fine electrified jet or droplets from a nozzle','DC or pulsed high voltage; precision stage','single-nozzle lab printers; multi-nozzle pilots'],
'inks: metal nanoparticles 10–60 wt%, polymers, sol–gels|conductivity 10⁻⁶–10⁻¹ S/m|viscosity 1–10 000 mPa·s','ambient printing; sintering 120–300 °C or photonic','jet much thinner than the nozzle|flow fL–pL per pulse|substrate charge build-up on insulators',
'0.2–3 kV|stand-off 10–500 µm|pulse 0.1–10 kHz for drop-on-demand','nozzle 0.1–30 µm inner diameter|stage accuracy sub-µm',
'voltage/pulse|stand-off|ink conductivity and viscosity|stage speed','line width|thickness|pitch|conductivity',[9,11,12,21,26],'high voltage|nanoparticle inks',
'Additive, mask-free and nearly waste-free; serial writing limits area.',[100, 5000],'writing speed 0.1–100 mm/s');
e(58,'N M 2 ox ce me sc pe hp L P I',[20,300],[5e-3,5e-2],[600,36000],
['Ar (Xe, Kr) with O₂ for reactive deposition','ion-beam-sputtered flux in high vacuum','gridded RF or Kaufman ion source','R&D tools to precision-optics production'],
'targets Ta, Si, Ti, Hf, metals 99.99 %+|O₂ near the substrate gives full oxidation','substrate 20–300 °C|low heat load','plasma confined to the source: low-pressure, collision-poor transport|sputtered atoms 5–20 eV|second (assist) source optional',
'beam energy 0.5–1.5 keV|beam current 0.1–1 A|neutraliser prevents charging','dense amorphous films, optical loss < 1–10 ppm|high compressive stress 0.1–0.6 GPa',
'beam energy|beam current|O₂ flow|assist energy','thickness (±0.1 %)|refractive index|absorption/scatter|stress',[7,10,12,16,21],'high voltage|vacuum hazards',
'Slow but the lowest-loss coatings available; long runs make uptime and target utilisation the levers.',null,'0.01–0.3 nm/s');
e(59,'N M 2 sc ox me ce hp L P I',[-100,600],[1e-5,1e-3],[1,3600],
['high vacuum beamline','energetic ions stopped in a solid','ion source, mass analyser and accelerator','R&D implanters to 300 mm production'],
'species B, P, As, Si, N, O, H, He, metals|dose 10¹¹–10¹⁷ ions/cm²|mass-analysed isotopic purity','wafer −100 to 600 °C (hot implants for SiC)|activation anneal 900–1100 °C afterwards|beam power heating','projected range 1 nm–several µm with straggle|channelling along crystal axes; 7° tilt|sputtering at high dose',
'energy 0.2 keV–3 MeV|beam current µA–30 mA|electrostatic or mechanical scanning','lattice damage and amorphisation above a critical dose|wafer charging neutralised by flood gun',
'ion species|energy|dose|tilt/twist','depth profile|peak concentration|damage level|buried-layer thickness',[11,16,17,21,23],'toxic source gases (AsH₃, PH₃, BF₃)|X-rays, high voltage',
'Precise and material-frugal; toxic gas use is cut by sub-atmospheric delivery sources.',null,'10–500 wafers/h by dose');
e(60,'M N 2 3 me sc ox ce ca L',[15,30],[1e-5,1e-3],[1,36000],
['high vacuum; optional precursor or etch-assist gas','focused ion beam on a solid surface','liquid-metal (Ga⁺) or gas-field (He⁺, Ne⁺), plasma Xe⁺ sources','single samples; fields of ~100 µm'],
'physical sputtering, gas-assisted etching (XeF₂, I₂) or deposition (Pt, W, C precursors)|Ga implantation ~10²⁰–10²¹ cm⁻³ near surface','ambient stage; cryo-stages for soft matter|local beam heating small in conductors','redeposition limits aspect ratio (≈ 5–10:1)|precursor gas flux via injection needle',
'Ga⁺ 5–30 keV|beam current 1 pA–65 nA|spot 5–500 nm|dwell 0.1–100 µs','amorphised damage layer 10–30 nm at 30 keV|drift and charging on insulators',
'beam current|ion dose|dwell/scan strategy|energy','feature width|depth|sidewall angle|damage',[5,12,16,21,23,29],'vacuum and high voltage|toxic precursor gases',
'Maskless and immediate for prototypes; far too slow for area production.',[5, 30],'milling ≈ 0.1–1 µm³/nC; serial');

// ── D. Lithography, cleanroom and precision nanofabrication ──
e(61,'M N 2 3 sc ox me ce pe hp L P I',[20,23],A,[3600,5e6],
['HEPA/ULPA-filtered air, ISO class 3–7, slight overpressure','controlled facility housing deposition, lithography and etch tools','air handling 50–600 air changes/h; facility utilities','university cleanrooms to gigafabs'],
'ultrapure water 18.2 MΩ·cm|airborne molecular contamination ppb|electronic-grade chemicals','air 20–23 °C ±0.1–1 K|humidity 40–45 % ±2–5 %','unidirectional downflow 0.3–0.5 m/s|particle limits by ISO 14644-1 class',
'electrostatic discharge control|yellow light (> 500 nm) in lithography bays','vibration criteria VC-D/E for lithography|gowning and material-transfer protocols',
'cleanliness class|temperature/humidity stability|vibration|process sequence','defect density|yield|overlay|line-width control',[5,12,21,23,26,27],'toxic, corrosive and pyrophoric process chemicals',
'Air handling dominates energy use; mini-environments, heat recovery and water reclaim give the largest savings. The time window spans single steps to full process flows.',null,'full flows take days–weeks');
e(62,'M N 2 sc me ox po pe L P I',[20,130],A,[30,3600],
['cleanroom air, yellow light','photoresist film exposed through a mask','Hg lamp (436/365 nm), KrF 248 nm or ArF 193 nm source','contact aligners to 300 mm scanners'],
'resist thickness 0.1–10 µm (to 100+ µm for SU-8)|aqueous developer TMAH 2.38 %|adhesion promoter HMDS','soft bake 90–115 °C, post-exposure bake 90–130 °C, 60–90 s','spin coat 1000–6000 rpm|developer puddle 30–90 s|immersion water layer for 193i',
'dose 10–500 mJ/cm²|resolution ≈ k₁λ/NA (NA up to 1.35)|depth of focus ≈ λ/NA²','substrate flatness and reflectivity (anti-reflective coatings)|overlay 2–500 nm by tool',
'wavelength/NA|dose|focus|resist and thickness','critical dimension|sidewall profile|overlay|pattern density',[5,12,21,23,26],'solvent vapour|TMAH is acutely toxic',
'Highly parallel; solvent, developer and ultrapure-water use dominate — track recycling and safer solvents help.',[40, 1000],'50–300 wafers/h');
e(63,'M N 2 0 1 me sc ox po pe hp L',[18,25],[1e-5,1e-3],[60,86400],
['high vacuum','electron-sensitive resist exposed by a focused beam','thermal field-emission gun, 10–100 keV','chips to 200 mm wafers and photomasks'],
'resists PMMA, ZEP, HSQ, CSAR at 30–500 nm|developer MIBK:IPA or TMAH|conductive layer on insulating substrates','stage and room ±0.1 K to limit drift','serial writing; write time ∝ dose × area / beam current|proximity effect from back-scattered electrons (range ~10–30 µm at 100 keV in Si)',
'10–100 keV|beam current 50 pA–100 nA|dose 100–1000 µC/cm² (PMMA at 100 keV), higher for HSQ|spot 2–10 nm','stitching error 5–20 nm between write fields|charging distorts patterns',
'dose|beam current/spot|acceleration voltage|resist and thickness','line width|pitch|edge roughness|placement accuracy',[5,12,21,23,25],'high voltage and X-rays|developer solvents',
'Maskless and flexible but energy- and time-intensive per area; reserve for masters, masks and prototypes.',[5, 20],'serial: mm²–cm² per hour for dense patterns');
e(64,'M N 2 3 po sc me ox pe lc L P I',[20,200],A,[10,600],
['air or evacuated imprint chamber (avoids trapped bubbles)','thermoplastic or UV-curable resist moulded by a stamp','heated press or UV lamp','chips to 300 mm step-and-repeat and roll-to-roll'],
'UV resists (acrylate, 1–20 mPa·s) or thermoplastics (PMMA)|anti-stick fluorosilane layer on stamp|residual layer removed by O₂ plasma','thermal NIL 50–100 K above T_g (100–200 °C), demould below T_g|UV-NIL at room temperature','squeeze flow fills cavities; time ∝ viscosity/pressure|residual layer 5–50 nm',
'UV dose 20–2000 mJ/cm² at 365–405 nm','imprint pressure 0.1–5 MPa (thermal), < 0.1 MPa (UV)|demoulding force and stamp wear|defectivity from particles',
'stamp pattern|pressure|temperature or UV dose|resist viscosity','feature size|aspect ratio|residual-layer uniformity|fidelity',[5,12,14,21,26],'UV exposure|resist sensitisers',
'Replicates one master thousands of times with little energy; the lowest-cost route to large-area sub-100 nm patterns.',[10, 50],'20–100 wafers/h; m/min roll-to-roll');
e(65,'M N 2 sc me pe hp I',[20,23],[1,10],[10,300],
['vacuum with ~5 Pa H₂ background','thin resist exposed by reflective 13.5 nm projection','laser-produced tin plasma source, ~250–500 W at intermediate focus','high-volume 300 mm fabs only'],
'chemically amplified or metal-oxide resists 20–40 nm thick|photon shot noise causes stochastic defects','optics and reticle thermally stabilised to mK','all-reflective Mo/Si multilayer optics (~70 % per mirror)|hydrogen keeps optics free of tin and carbon',
'wavelength 13.5 nm (92 eV)|NA 0.33 (≈ 13 nm half-pitch) or 0.55 (≈ 8 nm)|dose 30–80 mJ/cm²','pellicle and reticle defectivity|overlay ~1–2 nm',
'dose|NA/illumination|resist|focus','half-pitch|line-edge roughness|stochastic defect rate|overlay',[5,12,21,23,26],'high-power laser and hydrogen|tin debris',
'Roughly 1 MW electrical per scanner; replaces multi-patterning steps, which offsets part of the energy cost.',[8, 20],'> 100 wafers/h');
e(66,'M N 2 3 po sc me ox pe lc L P',[18,25],A,[1,600],
['air on a vibration-isolated table','photoresist exposed by interfering coherent beams','UV or visible laser (325, 355, 266, 193 nm)','cm² to m² gratings'],
'resist 50–1000 nm; anti-reflective coating suppresses standing waves|multi-exposure or multi-beam gives 2-D and 3-D lattices','stable to ±0.1 K during exposure','no mask; fringe stability better than period/10 needed|immersion shortens period',
'period = λ/(2 sin θ)|coherence length > path difference|dose 10–200 mJ/cm²','vibration and air-turbulence control|large-area uniformity limited by beam profile',
'wavelength|interference angle|dose|number of beams/exposures','period|duty cycle|lattice symmetry|uniform area',[5,12,15,21],'laser radiation',
'Maskless periodic patterning over large areas at low cost; limited to periodic designs.',[50, 500],'exposure seconds–minutes per field');
e(67,'M N 3 2 po co ce 3a pe L',[18,28],A,[60,172800],
['air, yellow light','liquid or solid photoresist polymerised at a laser focus','femtosecond laser 780 nm, ~100 fs, 80 MHz','structures µm³ to cm³; meso-scale with stitching'],
'acrylate, epoxy (SU-8) or hybrid (ORMOCER) resists with two-photon initiators|development in PGMEA/IPA|pyrolysis 900 °C gives glassy carbon (≈ 5× shrinkage)','ambient; local heating small below the damage threshold','voxel overlap and hatching/slicing 0.1–1 µm|oxygen quenching confines polymerisation|capillary collapse on drying — use critical-point drying',
'two-photon absorption ∝ intensity²|average power 5–50 mW at focus|NA 0.8–1.4|scan 0.1–100 mm/s','shrinkage 2–10 %|strut diameters 0.2–5 µm|adhesion to substrate',
'laser power|scan speed|objective NA|hatching/slicing distance','feature (voxel) size|lattice geometry|degree of conversion|stiffness',[5,12,14,15,21,25],'class 4 laser|resist sensitisers',
'Minimal material use and no masks; serial writing is slow, so reserve for unit cells and masters, then replicate.',[100, 300],'voxel rate 10³–10⁶ per second');
e(68,'M N 2 3 po sc me pe lc L P',[18,28],A,[60,86400],
['air, yellow light','photoresist or functional film scanned by a focused laser','diode/solid-state laser 355–405 nm (one-photon) or fs laser (multiphoton)','mm² to 200 mm and larger substrates'],
'positive or negative resists 0.5–100 µm|greyscale exposure for 2.5-D relief|direct ablation or photoreduction of functional films','ambient; thermal modes for phase-change or ablative writing','raster or vector scan; write time ∝ area/spot²',
'wavelength 355–405 nm|spot 0.3–5 µm|dose 10–500 mJ/cm²|intensity modulation for greyscale','autofocus on non-flat substrates|stage accuracy 10–100 nm',
'spot size/objective|dose|scan speed|greyscale levels','line width|relief depth|pattern area|sidewall profile',[5,12,14,15,21],'laser radiation|developer solvents',
'Maskless prototyping removes mask sets and their lead time.',[200, 1000],'1–1000 mm²/min by resolution');
e(69,'M N 2 sc ox me po pe L P I',[-20,150],[1,30],[10,3600],
['CF₄, CHF₃, SF₆, Cl₂, O₂, Ar','capacitively coupled plasma over a powered wafer electrode','13.56 MHz RF, 50–500 W','R&D chambers to production etchers'],
'fluorine or chlorine radicals form volatile products|fluorocarbon polymer passivates sidewalls|O₂ addition tunes selectivity','electrode cooled to 10–20 °C typical|wafer heating by ion bombardment','ion-assisted etching gives anisotropy|micro-loading and aspect-ratio-dependent etching|gas flow 10–200 sccm',
'DC self-bias 100–600 V|ion energy 100–500 eV|density 10⁹–10¹⁰ cm⁻³','mask erosion and faceting|ion damage and surface residues',
'RF power/bias|pressure|gas chemistry|time','etch depth|anisotropy/sidewall angle|selectivity|roughness',[5,7,11,12,21],'toxic and corrosive gases|RF',
'Dry pattern transfer avoids wet waste; fluorinated-gas abatement is essential.',null,'etch 10–500 nm/min; resolution set by the mask');
e(70,'M N 3 2 sc pe 3a L P I',[-110,40],[1,10],[60,14400],
['SF₆ etch alternating with C₄F₈ passivation (Bosch) or SF₆/O₂ (cryogenic)','high-density ICP plasma on silicon','ICP source 1–5 kW plus pulsed bias','100–300 mm wafers'],
'Bosch cycle: 1–10 s etch, 1–5 s passivation|selectivity to resist 50–100:1, to SiO₂ 150–300:1','chuck 0–40 °C (Bosch) or −110 °C (cryogenic)|He backside cooling','aspect-ratio-dependent etch rate|notching at buried oxide|loading effect with exposed area',
'ICP 1–5 kW at 13.56 MHz|low-frequency pulsed bias 10–100 W','scallops 20–200 nm|aspect ratio 20–50:1 (to ~100:1)|through-wafer etching needs carrier bonding',
'cycle times|ICP and bias power|pressure|temperature','depth|sidewall angle and scallop size|aspect ratio|uniformity',[5,14,15,21,23],'fluorinated greenhouse gases|RF power',
'SF₆ and C₄F₈ are potent greenhouse gases; high-efficiency abatement is the dominant environmental control.',[200, 2000],'etch 1–20 µm/min');
e(71,'N M 2 3 sc ox me pe lc L P I',[15,90],A,[5,7200],
['air, fume-extracted wet bench','liquid etchant (acid, base, or metal-assisted)','chemical; bath heating','beakers to automated wet benches'],
'KOH 20–40 wt% or TMAH 25 % for Si; (100):(111) selectivity 100–400:1|buffered HF for SiO₂ ~100 nm/min|metal-assisted etching: HF/H₂O₂ with Au or Ag','KOH 60–90 °C; BOE at 20–25 °C|rate roughly doubles per 10 K','reaction- or diffusion-limited; agitation matters|H₂ bubbles mask the surface|isotropic undercut ≈ etch depth',
'optional anodic bias for porous Si and alumina (1–200 mA/cm², 20–200 V)','stiction on drying of released structures|mask adhesion',
'etchant concentration|temperature|time|agitation','etch depth|undercut|pore diameter and depth|surface roughness',[6,7,9,11,21,27],'HF is acutely toxic|hot caustics',
'Cheap and parallel; chemical and rinse-water consumption and HF waste treatment are the burdens.',[500, 5000],'Si in KOH ≈ 1 µm/min at 80 °C; feature limit ≈ film thickness (far finer for anisotropic or metal-assisted etching)');
e(72,'M N 3 2 0 me ca ox 3a L',[15,30],[1e-4,1e-2],[1,36000],
['high vacuum with locally injected precursor gas','adsorbed precursor decomposed by a focused electron beam','SEM column 1–30 keV','single structures; µm-scale fields'],
'precursors MeCpPtMe₃, W(CO)₆, Co₂(CO)₈, TEOS|deposits often 15–60 at% metal in carbon matrix; purified by O₂/H₂O or annealing|etching with XeF₂ or H₂O','ambient stage; mild substrate heating changes adsorption','precursor-limited vs electron-limited regimes|surface diffusion replenishes adsorbates',
'1–30 keV|beam current 1 pA–10 nA|dwell 0.1–1000 µs|secondary electrons drive dissociation','free-standing 3-D wires and helices possible|drift limits long writes',
'beam current|dwell and refresh time|precursor flux|energy','feature width|height|metal content|3-D shape',[5,7,12,13,17,21],'toxic organometallics|X-rays',
'True 3-D nanoprinting without masks or resists; throughput suits prototypes and mask repair only.',[3, 30],'serial; ≈ 10⁻⁴–10⁻¹ µm³/nC, precursor-dependent');
e(73,'N M 2 0 po me sc ox bi L',[15,30],A,[1,86400],
['air at controlled humidity; UHV and cryogenic for atom manipulation','sharp probe tip acting on a surface','AFM/STM tip: mechanical, thermal, electrical or molecular-ink transfer','single tips to parallel cantilever arrays'],
'dip-pen inks (thiols, silanes, biomolecules)|local anodic oxidation of Si or Ti|thermally decomposable resists (PPA)','ambient; heated tips 300–1000 °C for thermal probe lithography','ink transport through the water meniscus; dot area ∝ dwell time|humidity 30–70 % controls meniscus',
'tip bias 5–20 V for oxidation|tunnelling current pA–nA','tip wear and drift|contact force nN–µN',
'dwell time/scan speed|humidity|tip bias or temperature|force','line width|dot size|pattern height|molecular coverage',[8,12,17,21,23],'low hazard',
'Almost no material or energy use; inherently slow unless massively parallel.',[10, 100],'writing 0.1–100 µm/s (to mm/s thermal)');
e(74,'M N 2 me po ce sc pe lc L P I',[15,30],A,[1,3600],
['air','surface directly ablated or melted by interfering pulsed beams','ns, ps or fs laser at 266–1064 nm','cm² samples to roll and sheet processing'],
'direct structuring of metals, polymers, ceramics — no resist|oxide and recast layers on metals','bulk near ambient; local melting/ablation at intensity maxima|thermal diffusion length sets minimum period for ns pulses (~1 µm in metals)','melt flow driven by recoil pressure and Marangoni forces|ablation plume removal',
'period = λ/(2 sin θ), 0.2–30 µm|fluence 0.1–5 J/cm²|pulse 10 fs–10 ns','structure depth up to ~period|surface hardening or heat-affected zone',
'fluence|period (angle)|pulse number|pulse duration','period|depth|wettability|optical diffraction and friction',[5,10,12,14,21,26],'class 4 laser|ablation fume',
'One-step, chemical-free surface texturing at industrial speed.',[200, 5000],'cm²/s-class; of order 1 m²/min reported');
e(75,'M N 2 3 sc me ox pe hp I P',[20,1050],[0.1,1.1e5],[3600,6e5],
['fab ambient with tool-specific vacuum and gases','sequences of deposition, etch and anneal using existing features as masks','combined ALD/CVD, plasma etch and RTP tools','high-volume semiconductor fabs'],
'spacer films SiO₂/SiN by ALD|self-aligned silicide: metal + Si reaction then selective strip|area-selective deposition chemistries','ALD spacers 50–400 °C|silicide anneals 250–850 °C|activation spikes to ~1050 °C','conformal deposition and anisotropic etch-back define spacers|pattern loading effects',
'plasma etch biases as in RIE/ICP','pitch halving (SADP) or quartering (SAQP)|pitch walking from mandrel CD error',
'mandrel CD|spacer thickness|etch anisotropy|anneal budget','final pitch and half-pitch|CD uniformity|line-edge roughness|overlay tolerance',[5,11,16,21,23,26],'process gases and wet chemicals of each step',
'Avoids extra lithography exposures at the cost of more deposition/etch steps; the window spans whole module flows.',[8, 40],'module flows take hours–days');
e(76,'M N 3 2 sc ox me ce hp L P I',[20,1100],[1e-3,1.1e5],[60,36000],
['cleanroom air, N₂ or vacuum','two polished solid surfaces brought into contact','press, heater, plasma activation; optional voltage','chips to 300 mm wafers'],
'hydrophilic direct bonding via –OH groups; plasma activation raises low-temperature strength|layer transfer by H⁺ implant ~5×10¹⁶ cm⁻² and splitting at 400–600 °C|adhesive, eutectic or metal thermocompression options','direct bond at room temperature then anneal 200–1100 °C|anodic 300–450 °C|eutectic Au–Sn ~280 °C, Cu–Cu 300–400 °C','trapped gas and particles cause voids|bond-wave propagation mm/s',
'anodic bonding 200–1000 V','surface roughness < 0.5 nm RMS and bow < tens of µm|thermal-expansion mismatch stress|alignment 0.1–2 µm',
'surface activation|anneal temperature|bonding pressure|alignment','bond energy (to ~2.5 J/m²)|void density|transferred-layer thickness|alignment accuracy',[8,10,14,21,22,23],'HF/piranha cleans|high voltage',
'Layer transfer reuses the donor wafer many times, saving high-value substrate material.',null,'5–30 wafer pairs/h');

// ── E. Mechanical, structural and additive manufacturing ──
e(77,'N 0 3 me ce ox co lc L P I',[-196,200],A,[1800,360000],
['sealed vial under Ar, air or liquid N₂ (cryomilling)','dry or wet powder charge with grinding media','planetary, attritor, shaker or stirred-media mills','g vials to tonne-scale attritors'],
'process-control agent 1–2 wt% (stearic acid, ethanol) limits cold welding|oxygen and nitrogen pick-up','vial bulk 20–200 °C; cryomilling at −196 °C|transient impact heating','repeated fracture and cold welding|wet milling adds slurry viscosity effects',
'','ball-to-powder ratio 5:1–20:1|200–600 rpm (planetary)|media wear contamination 0.1–several wt%',
'milling time|speed/energy|ball-to-powder ratio|media size','crystallite size 5–50 nm|particle size|lattice strain|contamination',[14,16,21,22,26],'pyrophoric fine metal powders|dust explosion',
'Solvent-free but energy-intensive with long times; stirred-media mills are several times more efficient than planetary mills.',null,'g/h (lab) to t/h (stirred media)');
e(78,'N 3 0 me co ce P I L',[20,300],A,[3600,360000],
['Ar or N₂ sealed','blended elemental or pre-alloyed powders under high-energy milling','attritor or planetary mill','lab to industrial ODS-alloy production'],
'elemental blends alloy by repeated welding and fracture|extended solid solubility, amorphisation|oxide dispersoids (Y₂O₃ 0.2–0.5 wt%)','bulk 20–300 °C|consolidation later by HIP, extrusion or SPS at 900–1200 °C','lamellar refinement to nm scale|steady state when welding balances fracture',
'','ball-to-powder ratio ~10:1|process-control agent 0–2 wt%|contamination from media and atmosphere',
'milling time|energy input|composition|process-control agent','grain size 10–50 nm|homogeneity|dispersoid size|phase',[6,14,16,21,22,26],'reactive powders|dust explosion',
'Makes alloys impossible by melting; milling energy is high, so maximise charge and use efficient attritors.',null,'batches kg–t over 10–100 h');
e(79,'N 3 me P I L',[20,500],A,[1,3600],
['air; lubricated tooling','bulk metal deformed under high hydrostatic pressure','hydraulic presses, rolling mills (ARB), forging','lab billets to continuous sheet and rod'],
'alloy purity and solutes set the steady-state grain size|precipitation during or after deformation','cold to warm working, below ~0.4 of melting temperature|adiabatic heating at high rates','no mass transport; dislocation cell formation then high-angle boundaries',
'','equivalent strain > 4–10|hydrostatic pressure 0.5–6 GPa suppresses cracking|accumulative roll bonding: 50 % reduction per cycle',
'accumulated strain|temperature|strain path|pressure','grain size 100–500 nm|yield strength|ductility|texture',[14,16,18,21,22],'high-force press hazards',
'Strengthens by grain refinement without alloying additions — less critical raw material — at significant forming energy.',null,'continuous variants m/min');
e(80,'N 3 me P L I',[20,450],A,[10,3600],
['air; MoS₂ or graphite lubricant','billet pressed through an angled channel','hydraulic press','cm billets to continuous ECAP-Conform'],
'Al, Cu, Mg, Ti and steels; hard-to-deform alloys pressed warm','room temperature (Al, Cu) to 200–450 °C (Mg, Ti)|die heating','simple shear at the channel intersection; cross-section unchanged',
'','channel angle 90–120°; strain ≈ 1 per pass at 90°|4–8 passes; routes A, B_C, C|punch pressure 0.2–1.5 GPa; back-pressure aids ductility|ram speed 1–10 mm/s',
'number of passes|route|die angle|temperature','grain size 0.2–1 µm|strength|texture|homogeneity',[14,16,18,21,22],'high-pressure tooling failure',
'Bulk nanostructuring with no material loss apart from billet ends; die life and lubrication are the cost drivers.',null,'minutes per pass');
e(81,'N 3 2 me co L',[-196,400],A,[60,3600],
['air','thin disc compressed and twisted between anvils','press with rotating anvil','discs 10–40 mm diameter, ~1 mm thick; lab scale'],
'pure metals, alloys, immiscible systems and powders|strain-induced phase transformations and supersaturation','usually ambient; cryogenic or up to 400 °C possible|temperature rise a few tens of K at 1 rpm','no bulk transport; forced mechanical mixing at very high strain',
'','pressure 1–10 GPa (typically 2–6)|shear strain γ = 2πNr/h; 5–10 turns at ~1 rpm|strain gradient from centre to edge',
'pressure|number of turns|rotation speed|temperature','grain size 10–200 nm|hardness|phase|homogeneity radius',[14,16,17,21,22],'high stored elastic energy in tooling',
'The finest bulk grain sizes of any deformation route, but sample size confines it to research and small parts.',null,'minutes per disc');
e(82,'N 2 0 1 ca ox co lc L P I',[0,60],A,[60,36000],
['air; ice bath or jacket cooling','liquid dispersion under ultrasonic cavitation','bath (40 kHz) or probe (20 kHz) sonicator','vials to continuous flow cells'],
'solvent surface tension ≈ 40 mJ/m² (NMP, DMF) or aqueous surfactant for 2-D materials|solids 0.1–50 g/L','bulk kept below ~40 °C; solvent degrades with heating','cavitation shear exfoliates and de-agglomerates|scission shortens nanotubes and flakes with time',
'20–40 kHz|probe 100–750 W; energy density 0.1–10 kJ/mL','probe-tip erosion|centrifugation selects size afterwards',
'energy density|time|solvent/surfactant|concentration','flake thickness and lateral size|dispersion concentration|agglomerate size|defect density',[8,9,14,21,26],'aerosols and noise|toxic solvents (NMP)',
'Simple but inefficient (yields often < 5 %); prefer green solvents and high-shear mixing at scale.',null,'exfoliation yields 0.01–1 g/h');
e(83,'M N 3 me co 3a pe P I L',[600,3000],[1e5,1.1e5],[3600,6e5],
['Ar or N₂, O₂ < 0.1 %','metal powder bed selectively melted by laser','fibre laser 1070 nm, 100–1000 W, 1–12 beams','build volumes 100 mm cube to ~1 m'],
'gas-atomised powder 15–45 µm (Ti-6Al-4V, 316L, AlSi10Mg, IN718)|oxygen and moisture pick-up|vaporisation loss of volatile elements','melt pool above the alloy liquidus (≈ 600 °C Al alloys to > 2000 °C peak); plate preheat 80–500 °C|cooling 10⁵–10⁶ K/s|residual stress up to yield','Marangoni convection and keyhole regime|shield-gas cross-flow removes spatter|powder recoating',
'spot 50–100 µm|scan 0.3–2 m/s|absorptivity 0.3–0.7','layer 20–60 µm, hatch 80–120 µm|energy density 40–120 J/mm³|supports below ~45° overhang|minimum strut ≈ 0.15–0.3 mm',
'laser power|scan speed|hatch and layer thickness|scan strategy','density (> 99.5 %)|strut diameter|surface roughness Ra 5–15 µm|microstructure',[5,10,14,16,21,22],'reactive metal powder fire/explosion|inert-gas asphyxiation',
'Near-net shape with high powder reuse; specific energy is high (tens of kWh/kg), justified by lightweight lattices.',[100000, 300000],'5–100 cm³/h');
e(84,'M N 3 me 3a pe P I L',[600,3000],[1e-3,1],[3600,6e5],
['vacuum with ~0.2 Pa He','pre-sintered hot metal powder bed melted by electron beam','electron gun 60 kV, 3–6 kW','build volumes up to ~350 mm diameter'],
'powder 45–105 µm (Ti-6Al-4V, TiAl, CoCr, Cu, Ni superalloys)|vacuum protects reactive alloys|Al evaporation loss','bed held at 600–1100 °C, giving low residual stress; melt pool above liquidus|slow cool of the whole build','powder lightly sintered to avoid charge-driven "smoke"|deep, fast-scanned melt pool',
'beam current 1–50 mA|magnetic deflection up to km/s|multi-spot melting','layer 50–100 µm|rougher surface Ra 20–35 µm|minimum strut ≈ 0.3–1 mm|sintered powder removal from fine lattices',
'beam current|scan speed|bed temperature|focus offset','density|strut diameter|microstructure (α/β, γ-TiAl)|roughness',[5,10,14,16,21,22],'X-rays|reactive powder',
'Hot-bed process needs little stress relief and reuses powder well; vacuum and preheat add energy.',[300000, 1000000],'50–100 cm³/h');
e(85,'M N 3 po co 3a pe P I L',[150,380],A,[3600,2e5],
['N₂ with a few % O₂ or less','polymer powder bed sintered by laser','CO₂ laser 10.6 µm, 30–100 W (or diode/fibre with absorber)','build volumes 200–700 mm'],
'PA12, PA11, TPU, PEEK powders 20–80 µm|refresh rate 30–50 % virgin powder|filled grades (glass, carbon, Al)','bed held just below melting: 170–180 °C for PA12, up to ~340–380 °C for PEEK|slow cool-down 10–40 h limits warpage','viscous sintering and coalescence|powder supports the part — no support structures',
'spot 0.3–0.5 mm|scan 2–10 m/s','layer 100–150 µm|minimum wall 0.5–0.8 mm|porosity 2–5 %, anisotropy in z',
'laser energy density|bed temperature|layer thickness|powder refresh ratio','density|strut/wall size|elongation at break|surface roughness',[5,10,14,15,21,26],'combustible powder dust|hot surfaces',
'Support-free lattices; aged powder is the main waste stream — maximise nesting density and refresh efficiency.',[400000, 800000],'10–40 mm build height per hour');
e(86,'M N 3 po co ce 3a pe lc L P I',[20,40],A,[1800,172800],
['air','liquid photopolymer vat cured by a scanned UV laser','UV laser 355 or 405 nm, 0.1–2 W','desktop to 1.5 m industrial vats'],
'acrylate/epoxy resins with photoinitiator 0.5–5 wt%|ceramic- or nanoparticle-filled resins to 40–60 vol%|wash in IPA; UV post-cure','vat 20–40 °C; post-cure 60–80 °C','recoating of viscous resin (0.1–5 Pa·s)|oxygen inhibition at free surfaces|drain holes for trapped resin',
'cure depth C_d = D_p ln(E/E_c)|spot 80–140 µm','layer 25–100 µm|cure shrinkage 2–8 % and warpage|supports required',
'exposure energy|layer thickness|resin (D_p, E_c)|orientation','feature size|green strength/modulus|surface finish|dimensional accuracy',[5,12,14,15,21,22],'skin-sensitising uncured resin|solvent wash',
'Fine architected polymers at low energy; resin and wash solvent are hazardous waste — cure before disposal.',[50000, 200000],'10–30 mm build height per hour');
e(87,'M N 3 po co ce 3a pe lc L P I',[20,45],A,[600,86400],
['air; oxygen-permeable window for continuous variants','photopolymer cured layer-wise by projected images','UV LED 385/405 nm with DMD projector','micro-scale (PµSL) to 300+ mm build areas'],
'low-viscosity resins (< 1–3 Pa·s) with absorber dye for z-resolution|filled, elastomeric and multi-material resins','vat 20–45 °C; exothermic cure in thick sections','resin reflow under the part limits speed|oxygen dead zone 20–100 µm enables continuous printing',
'irradiance 1–20 mW/cm²|pixel 25–100 µm (2–10 µm with projection micro-stereolithography)','layer 10–100 µm|peel forces on fine struts|post-cure shrinkage',
'exposure time/irradiance|layer thickness|pixel size|absorber content','feature size|strut diameter|conversion/modulus|accuracy',[5,12,14,15,21,25],'sensitising resin|UV light',
'Whole layers cure at once, so energy per part is very low; manage resin and solvent waste.',[2000, 100000],'10–100 mm/h (continuous: several hundred mm/h)');
e(88,'M 3 po co 3a pe lc L P I',[180,420],A,[600,6e5],
['air; heated chamber for engineering polymers','thermoplastic filament extruded through a hot nozzle','resistive hot end; motion gantry','desktop to room-sized printers'],
'PLA, ABS, PETG, TPU, PA, PEEK; fibre- or particle-filled grades|moisture < 0.02–0.05 % for hygroscopic polymers','nozzle 180–300 °C (PEEK 360–420 °C)|bed 20–120 °C; chamber to 90–200 °C|interlayer weld strength set by interface temperature history','shear-thinning melt flow; die swell|bridging and overhang limits',
'','nozzle 0.2–0.8 mm; layer 0.05–0.4 mm|anisotropy: z-strength 40–80 % of in-plane|warpage from crystallisation shrinkage',
'nozzle temperature|layer height|print speed|infill pattern/density','strut/wall size|porosity|stiffness and anisotropy|surface finish',[5,10,14,15,21,28],'ultrafine particle and VOC emissions|hot parts',
'Lowest-cost route to macroscale auxetic and cellular designs; recycled and bio-based filaments are available.',[200000, 800000],'5–100 cm³/h');
e(89,'M N 3 co ce po me ca 3a pe lc L P',[15,40],A,[60,86400],
['air; humidity or oil bath controls drying','yield-stress ink extruded as a continuous filament','pneumatic, screw or piston extrusion','micro-nozzles to large-format gantries'],
'colloidal, sol–gel, polymer, graphene oxide, metal or ceramic inks at 20–60 vol% solids|rheology modifiers','ambient printing; then drying, curing or sintering (to 1600 °C for ceramics)','shear-thinning: viscosity 10–10⁴ Pa·s at rest|storage modulus 10³–10⁵ Pa and yield stress 10–1000 Pa hold spanning filaments',
'','extrusion pressure 10–700 kPa|spans up to ~10–20 filament diameters|drying shrinkage and cracking',
'ink rheology|nozzle diameter|pressure/flow|print speed','filament diameter|lattice spacing|density after sintering|filler alignment',[5,8,9,14,21,25],'solvent vapour|nanoparticle inks',
'Room-temperature printing of almost any material with little waste; drying and sintering dominate energy.',[1000, 400000],'print speed 1–50 mm/s');
e(90,'M N 2 me po co ca cf L P',[15,80],A,[60,36000],
['air or N₂ sheath gas','focused aerosol stream of ink droplets','ultrasonic or pneumatic atomiser; sheath-gas focusing','single-nozzle R&D to multi-nozzle production'],
'nanoparticle inks (Ag, Au, Cu), polymers, dielectrics, CNT|viscosity 1–1000 mPa·s|particle size < 0.5 µm','platen 20–80 °C; sintering 120–250 °C, laser or photonic','droplets 1–5 µm|aerosol 10–50 sccm, sheath 20–100 sccm|overspray at the line edge',
'','stand-off 1–5 mm prints over steps and 3-D surfaces|line thickness 0.1–5 µm per pass',
'sheath-to-aerosol flow ratio|stage speed|nozzle size|ink loading','line width 10–100 µm|thickness|resistivity|edge definition',[9,11,12,21,26],'nanoparticle aerosol|solvent vapour',
'Additive and conformal with minimal ink waste; silver use dominates cost — copper inks reduce it.',[10000, 100000],'writing 1–200 mm/s');
e(91,'M N 2 me po ox co sc lc L P I',[15,70],A,[10,36000],
['air','picolitre droplets jetted onto a substrate','piezoelectric or thermal drop-on-demand heads','single nozzles to industrial multi-head lines'],
'viscosity 2–20 mPa·s, surface tension 25–35 mN/m|printability Z = 1/Oh between 1 and 10|particles < 1/50 of nozzle diameter','head 20–70 °C; substrate heated to control drying|sintering 120–250 °C or photonic','drop spreading by Weber/Reynolds numbers|coffee-ring effect; Marangoni suppression with co-solvents',
'piezo waveform 10–100 V at 1–50 kHz','drop volume 1–80 pL|substrate surface energy patterning sharpens lines|nozzle clogging',
'drop spacing|waveform|ink formulation|substrate temperature and wetting','line width 20–100 µm|film thickness|conductivity|pattern fidelity',[8,9,11,12,21,26],'solvent vapour|nanoparticle inks',
'Digital, mask-free and material-efficient; scales to roll-to-roll.',[20000, 100000],'up to m²/min industrially');
e(92,'M N 3 po co bi 3a lc L P',[4,250],A,[60,86400],
['air; sterile enclosure for biomaterials','paste, melt or gel extruded through a fine nozzle','pneumatic, piston or screw-driven dispenser','benchtop printers; multi-head systems'],
'polymer melts, pastes, hydrogels, filled composites|crosslinking by cooling, ions (Ca²⁺), UV or solvent loss','hydrogels 4–37 °C; thermoplastic melts 60–250 °C','viscosity 0.03–10⁴ Pa·s|shear rate 10–10³ s⁻¹ in the nozzle|filament merging and sagging',
'','nozzle 50–1000 µm|pressure 10–700 kPa|strand pitch defines pore size',
'pressure/flow|nozzle diameter|temperature|print speed','strand diameter|pore size|porosity|stiffness',[5,9,14,21,22],'hot melts|sensitising monomers',
'Low-energy, low-waste multimaterial structuring.',[50000, 1000000],'print speed 1–50 mm/s');
e(93,'M N 2 po me ox co cf pe lc P I',[20,200],[0.1,1.1e5],[0.1,600],
['air, N₂, or vacuum web coaters','continuous flexible web through coating, imprinting, curing and drying stations','web drives, dryers, UV/IR lamps, vacuum deposition zones','pilot lines 0.3 m wide to production > 1.5 m'],
'functional inks, UV-imprint resins, sputtered or evaporated layers|solvent drying load','web limited by substrate: PET ~150 °C, PEN ~180 °C, polyimide > 300 °C|drying ovens 60–150 °C','slot-die, gravure or flexographic coating; wet thickness 1–100 µm|drying residence seconds–minutes',
'UV cure 100–2000 mJ/cm²|corona/plasma pre-treatment','web tension 50–500 N/m|registration ±5–50 µm|roller defects repeat periodically',
'web speed|tension|coating gap/flow|cure dose','layer thickness|pattern fidelity|registration|defect density',[5,9,10,12,21,26],'solvent vapour (explosion limits)|nip points',
'Lowest cost and energy per m² for nanostructured films; solvent recovery or water-borne inks are the key measures.',[50, 10000],'0.1–100 m/min');
e(94,'M 3 2 co me po ce sc 3a pe L P I',[15,30],A,[1,86400],
['air or cleanroom','discrete parts positioned and joined','robots, micro-grippers, pick-and-place, transfer printing','single devices to surface-mount production lines'],
'adhesives, solders (reflow to ~260 °C), or mechanical interlock|surface cleanliness for bonding','ambient handling; local heating for bonding','no bulk transport; van der Waals and capillary adhesion dominate below ~100 µm part size',
'machine-vision alignment|electrostatic or vacuum gripping','placement accuracy 1–25 µm (sub-µm for transfer printing)|joint compliance alters resonator behaviour|tolerance stack-up across many cells',
'placement accuracy|joining method|part tolerance|assembly sequence','unit-cell position|periodicity error|joint stiffness|yield',[5,14,15,19,20,26],'robot pinch points|solder fume',
'Combines dissimilar materials without wasteful subtractive steps; design for disassembly aids reuse.',[1000, 50000],'1 to > 10⁴ placements per hour');
e(95,'M 3 2 po me co 3a pe lc L P',[15,150],A,[60,86400],
['air','sheet material cut, creased and folded','laser, blade or die cutting; manual, robotic or stimulus-driven folding','mm-scale films to architectural panels'],
'paper, polymer films, metal foils, composites, shape-memory polymers|active layers (hydrogel, bimorph) for self-folding','ambient; self-folding triggered at 60–150 °C or by light/moisture','not applicable (solid sheets); humidity affects paper and hydrogels',
'optional electrostatic, magnetic or light actuation','sheet thickness 10 µm–mm|crease stiffness vs panel bending stiffness|cut-pattern stress concentrations and fatigue at hinges',
'cut/crease pattern|sheet thickness|hinge stiffness|fold angle','deployed shape|Poisson\'s ratio|stretchability|multistability',[5, 14, 19, 20, 21, 22],'laser cutting fume|sharp edges',
'Flat-pack fabrication minimises material and transport volume; recyclable mono-material sheets are feasible.',[10000, 1000000],'cutting 10–1000 mm/s');
e(96,'M 3 2 me po ce co pe 3a L P I',[15,30],A,[60,86400],
['air with coolant or dielectric fluid','solid workpiece shaped by micro-milling, micro-EDM, laser or waterjet','high-speed spindles, EDM generators, ultrafast lasers','job-shop to production machining centres'],
'metals, polymers, ceramics, silicon, glass|coolant/dielectric chemistry|recast and heat-affected layers (EDM, ns-laser)','shop held at 20 ±0.5–1 K for accuracy|ultrafast lasers minimise heat-affected zone','chip evacuation and flushing of the cut|debris redeposition',
'micro-EDM pulses 0.1–10 µJ|fs–ps laser fluence 0.1–10 J/cm²','tools 25–1000 µm; spindle 20 000–200 000 rpm|tolerance ±1–5 µm; Ra 0.05–0.4 µm|burrs and tool deflection limit thin walls',
'tool/beam size|feed and speed|depth of cut|fixturing','feature size|tolerance|surface roughness|wall aspect ratio',[5,14,15,19,21,23],'swarf, coolant mist|laser radiation',
'Subtractive, so material yield is low for open lattices; recycle swarf and use minimum-quantity lubrication.',[1000, 100000],'removal 0.01–100 mm³/min');

// ── F. Biological and bioinspired synthesis ──
e(97,'N 0 me ox sc bi lc L P',[20,40],A,[14400,432000],
['aerobic or anaerobic, sterile','aqueous culture medium with living bacteria','metabolic; shaking/stirring and aeration','shake flasks to m³ fermenters'],
'metal salt 0.1–5 mM (toxicity-limited)|pH 5–9|carbon source and nutrients|reduction by reductases or secreted metabolites','25–37 °C ±0.5 K (mesophiles); thermophiles to 60–70 °C','oxygen transfer k_La 10–200 h⁻¹|stirring 100–500 rpm|intracellular products need cell lysis',
'','sterility|downstream separation from biomass',
'metal concentration|pH|temperature|strain','particle size 2–100 nm|shape|biogenic capping layer|yield',[6,7,8,21,27,28],'biosafety of strains|metal-laden biomass',
'Ambient, aqueous and reagent-light; slow, with broader size distributions and biomass waste to manage.',null,'mg–g per litre per day');
e(98,'N 0 me ox co bi lc L P',[20,35],A,[86400,604800],
['aerobic, sterile','aqueous broth or cell-free fungal filtrate','metabolic; aeration and agitation','flasks to industrial fermenters'],
'metal salt 0.5–5 mM|pH 4–8|extracellular enzymes and proteins reduce and cap (e.g. Fusarium, Aspergillus, Trichoderma)','25–30 °C typical','high-viscosity mycelial broths limit oxygen transfer|cell-free filtrate route simplifies separation',
'','sterility|protein corona on products',
'filtrate/biomass ratio|metal concentration|pH|temperature','particle size 5–100 nm|dispersity|capping proteins|yield',[6,7,8,21,27,28],'spores and allergens|mycotoxin-producing strains',
'High secretion makes extracellular synthesis and separation easier than with bacteria; dispose of biomass safely.',null,'mg–g per litre per day');
e(99,'N 0 me ox bi lc L P',[15,35],A,[86400,1.2e6],
['air enriched with 0.04–5 % CO₂','illuminated aqueous culture of micro- or macroalgae (or extracts)','light (photosynthesis); bubbling or paddle mixing','flasks to photobioreactors and open ponds'],
'metal salt 0.1–2 mM|pH 6–9|nutrients N, P|pigments, polysaccharides and proteins reduce and stabilise','20–30 °C','light attenuation limits depth to cm|gas exchange and mixing avoid settling',
'irradiance 50–300 µmol photons m⁻² s⁻¹|photoperiod e.g. 12 h:12 h','biofouling of reactor walls|harvesting by centrifugation or flocculation',
'metal concentration|light intensity|pH|biomass density','particle size|shape|biogenic coating|yield',[6,7,8,21,27,28],'low; metal-laden biomass disposal',
'Solar-driven and CO₂-fixing, and can couple to metal recovery from wastewater; productivity is low.',null,'mg per litre per day');
e(100,'N 0 me ox bi lc L P',[20,95],A,[60,86400],
['air','aqueous plant extract mixed with metal salt','stirring; mild heating','beakers to hundreds of litres'],
'metal salt 0.5–10 mM (AgNO₃, HAuCl₄, Zn, Cu, Fe salts)|polyphenols, flavonoids, sugars and terpenoids act as reductant and capping agent|pH 4–11','25–95 °C; oxides often need calcination 300–500 °C','fast mixing gives narrower sizes|extract composition varies with source and season',
'','batch-to-batch variability of extract|purification by centrifugation',
'extract:salt ratio|pH|temperature|time','particle size 5–100 nm|shape|surface chemistry|stability',[6,7,8,21,27,28],'low; silver and copper ecotoxicity',
'Renewable reductants and water solvent; standardise extracts to control reproducibility.',null,'g per batch in minutes–hours');
e(101,'N 0 2 ox me sc bi L',[20,60],A,[600,86400],
['air or controlled gas','buffered aqueous enzyme solution','biocatalysis; gentle mixing','µL assays to litre reactors'],
'enzyme nM–µM (urease, silicatein, glucose oxidase, reductases, peroxidase)|substrate and cofactor (NADH/NADPH) supply|buffer pH 5–9','at or below enzyme optimum (25–60 °C); denaturation above','Michaelis–Menten kinetics set supersaturation release rate|immobilised enzymes localise mineral formation',
'','enzyme stability and reuse|product inhibition',
'enzyme concentration|substrate concentration|pH|temperature','particle size|crystal phase|spatial localisation|morphology',[6,7,8,17,21],'low hazard',
'Mild and selective; enzyme cost is the main barrier, so immobilise and reuse.',null,'turnover-limited; mg per hour');
e(102,'N 0 1 ox me sc bi hp L',[4,40],A,[600,172800],
['air or N₂ for Fe²⁺ chemistry','buffered protein solution; cages, fibres or S-layers as templates','chemical; gentle stirring','mL to litres'],
'protein µM (ferritin, virus capsids, collagen, silk)|metal ions dosed in steps (e.g. Fe²⁺ + oxidant)|pH 6–9','4–40 °C, below denaturation','ion entry through protein channels|nucleation on charged residues inside the cage',
'','cage cavity defines size (ferritin 8 nm core, 12 nm shell)|protein integrity',
'metal:protein ratio|pH|addition rate|temperature','core size (monodisperse)|composition|crystal phase|array ordering',[3,7,8,17,21],'low hazard',
'Uniform particles under benign conditions; protein supply limits scale.',null,'mg per batch');
e(103,'N M 2 3 0 bi me pe hp L',[4,95],A,[3600,259200],
['air, sealed tubes','buffered aqueous solution of DNA scaffold and staples','thermal cycler annealing ramp','µL tubes; mg–g with biotechnological strand production'],
'scaffold 5–50 nM (M13, 7249 nt) with 5–10× staple excess|Mg²⁺ 10–20 mM in Tris–acetate–EDTA|nanoparticle attachment via thiolated DNA','denature at 65–95 °C then cool to 20–25 °C over 1–72 h|storage at 4 °C','diffusion-limited hybridisation; purification by filtration, gel or PEG precipitation',
'optional electric-field or lithographic placement on surfaces','structures 10–100 nm; denature above ~55–60 °C or at low Mg²⁺|nuclease sensitivity',
'sequence design|Mg²⁺ concentration|annealing ramp|staple excess','shape|addressable site positions (≈ 6 nm pitch)|yield|particle spacing',[3,4,5,12,17,25],'low hazard',
'Programmable with atom-economical chemistry; DNA cost is falling but still limits large volumes.',[6, 10],'folding 1–72 h; nM-scale yields');
e(104,'N 1 3 bi po lc L P',[4,80],A,[60,604800],
['air','aqueous or mixed-solvent peptide solution','self-assembly triggered by pH, salt, temperature or solvent switch','µL to litres'],
'peptide 0.1–20 mg/mL above the critical aggregation concentration|pH 2–11 and ionic strength control charge|β-sheet, α-helical coiled-coil or amphiphile designs','ambient; annealing 60–80 °C then cooling improves order','nucleation–elongation kinetics with lag phase|gelation above ~0.5–2 wt%',
'','hydrogel storage modulus 0.1–100 kPa|shear alignment of fibres',
'peptide sequence|concentration|pH|ionic strength','fibre diameter 5–20 nm|length|gel stiffness|morphology',[8,14,17,21,27],'low hazard',
'Biodegradable, water-processed building blocks; solid-phase peptide synthesis is solvent-heavy, so prefer short sequences.',null,'assembly seconds–days');
e(105,'N 0 3 bi po L P I',[20,65],A,[60,14400],
['air or N₂; sterile for pharmaceutical use','aqueous buffer with lipid in ethanol or as a dry film','hydration, extrusion, sonication or microfluidic mixing','mL to continuous GMP lines'],
'lipid 1–50 mM (phospholipid, cholesterol, PEG-lipid, ionisable lipid)|aqueous:ethanol flow-rate ratio ~3:1|pH 4 for nucleic-acid loading then buffer exchange','above the lipid transition temperature (DPPC 41 °C)|storage 2–8 °C or frozen','mixing time < 10 ms sets nanoparticle size|total flow 1–20 mL/min per mixer|tangential-flow filtration removes ethanol',
'','extrusion through 50–200 nm membranes, 10–21 passes, up to several MPa|high-pressure homogenisation 50–150 MPa',
'lipid composition|flow-rate ratio|total flow rate|temperature','diameter 50–200 nm|polydispersity|lamellarity|encapsulation efficiency',[8,9,17,21,26,27],'ethanol flammability',
'Aqueous and mild; ethanol recovery and cold-chain energy are the main footprint items.',null,'mL/min per mixer; L/h in parallel');
e(106,'N M 3 2 ce ox co bi lc L P',[4,40],A,[3600,2.4e6],
['air; CO₂ or NH₃ vapour diffusion for carbonates','aqueous solution with organic additives or matrices','chemical; slow diffusion or titration','vials to coating baths'],
'Ca²⁺ 1–20 mM with carbonate or phosphate|simulated body fluid pH 7.4|acidic polymers (polyaspartate 10–100 µg/mL) stabilise amorphous precursors','20–40 °C (37 °C for body-fluid routes)','slow supersaturation rise through gas or gel diffusion|infiltration of liquid-like precursor into collagen or templates',
'','organic matrix controls orientation and toughness|layered brick-and-mortar structures',
'supersaturation|additive type and concentration|pH|template','crystal polymorph|orientation|hierarchical structure|mineral content',[6,7,8,14,16,21],'low hazard',
'Ambient-temperature ceramics from abundant ions; growth is slow.',null,'µm per day');
e(107,'N M 3 1 po bi co 3a L P',[-80,100],[1,1.1e5],[3600,259200],
['air or sterile; vacuum during freeze-drying','polymer solution, melt or hydrogel shaped into a porous scaffold','electrospinning, freeze-drying, phase separation, porogen leaching, 3-D printing','mm implants to sheet production'],
'PCL, PLA, PLGA, collagen, gelatin, chitosan, bioceramics|porogen (salt, sugar) 70–90 vol%|crosslinkers (genipin, EDC)','freezing −20 to −80 °C; primary drying below ~100 Pa|melt processing 60–100 °C (PCL)','ice-crystal templating sets pores|residual solvent removal|interconnectivity governs nutrient transport',
'','porosity 70–95 %; pores 100–400 µm for bone|compressive modulus kPa–hundreds of MPa|sterilisation changes properties',
'polymer concentration|freezing rate|porogen size|fibre diameter','porosity|pore size and interconnectivity|stiffness|degradation rate',[5,9,14,21,22,27],'organic solvents (chloroform, HFIP)',
'Use biodegradable polymers with benign solvents (water, acetic acid); freeze-drying is the energy hotspot.',null,'batches in hours–days');
e(108,'N M 3 bi po co 3a L P',[4,37],A,[60,14400],
['sterile air; culture at 5 % CO₂ and 95 % humidity','cell-laden hydrogel bioink','extrusion, inkjet, laser-assisted or light-based bioprinting','mm–cm constructs'],
'bioinks: alginate 1–4 wt%, GelMA 5–15 wt%, collagen, fibrin|cell density 10⁶–10⁸ per mL|ionic (CaCl₂) or photo-crosslinking','cartridge 4–37 °C; bed 4–37 °C','shear-thinning bioink; viability falls as nozzle shear stress rises (kPa range)|diffusion limit ~200 µm without perfusion',
'photo-crosslinking 365–405 nm at 5–20 mW/cm²','nozzle 100–400 µm; pressure 5–200 kPa|construct stiffness 0.1–100 kPa',
'bioink rheology|pressure/nozzle|print speed|crosslinking dose','strand width|cell viability (> 80–90 %)|construct stiffness|pore architecture',[5,9,14,21,22,27],'biosafety of cells|photoinitiator cytotoxicity',
'Very small material volumes; single-use sterile consumables dominate waste.',[50000, 500000],'print speed 1–30 mm/s');

// ── G. Self-assembly, interfacial and soft-matter environments ──
e(109,'M N 3 2 po ox co pe lc L P',[15,65],A,[3600,6e5],
['air at controlled humidity and temperature','monodisperse colloidal suspension','evaporation, sedimentation, shear or confinement','cm² films to wafer scale'],
'silica or polymer spheres 100–1000 nm, polydispersity < 5 %|volume fraction 0.1–5 % (vertical deposition)|ionic strength < 1 mM for charge-stabilised crystals','20–65 °C, ±0.5 K','convective assembly at the meniscus|Péclet number: sedimentation vs Brownian motion|capillary forces lock the lattice',
'','drying cracks from shrinkage|infiltration and inversion give inverse opals',
'particle size and dispersity|concentration|evaporation rate|temperature','lattice constant|number of layers|domain size|stop-band wavelength',[5,8,9,12,15,21],'low hazard',
'Bottom-up photonic structures with almost no energy input; defects and cracks limit device use.',null,'film growth 0.1–10 µm/min');
e(110,'N M 2 3 ox po co pe lc L P I',[15,80],A,[10,86400],
['air at controlled humidity','drying sol or nanoparticle dispersion film','solvent evaporation during dip, spin, spray or blade coating','cm² to roll-to-roll coating'],
'surfactant or block-copolymer template (CTAB, Pluronic P123/F127) passes the critical micelle concentration on drying|silica or titania sol|nanocrystals with ligands in toluene or hexane','ambient coating; template removal 350–500 °C or by extraction','evaporation front concentrates the solute|Marangoni and coffee-ring flows|withdrawal 0.1–10 mm/s sets thickness',
'','relative humidity 20–80 % selects the mesophase|shrinkage on calcination 20–40 %',
'evaporation rate|humidity|template:precursor ratio|withdrawal speed','pore size 2–30 nm|mesophase symmetry|film thickness|superlattice order',[8,9,10,17,21],'flammable solvents',
'Fast and coating-line compatible; recover solvent and prefer extractable templates over burn-out.',null,'ordering in seconds–minutes');
e(111,'N M 2 po me ca bi hp L',[10,30],A,[1800,28800],
['clean air, vibration-isolated trough','monolayer at the air–water interface transferred to a solid','moving barriers and dipper','cm² substrates; continuous roll variants exist'],
'amphiphiles or ligand-capped nanoparticles spread from chloroform or hexane|ultrapure water subphase; pH and ions tune packing','subphase 15–25 °C ±0.1 K','compression isotherm: gas, liquid and solid phases|transfer ratio ≈ 1 for good films',
'','surface pressure 10–40 mN/m at transfer|barrier speed 1–10 mm/min; dipping 1–10 mm/min|collapse above critical pressure',
'surface pressure|dipping speed|subphase composition|number of layers','packing density|layer number|orientation|film thickness (monolayer precision)',[8,12,17,21,23],'chlorinated spreading solvents',
'Monolayer control with negligible material use; slow and small-area.',null,'1–10 mm/min per layer');
e(112,'M N 2 3 po sc ox pe lc L P I',[20,250],[1e2,1.1e5],[60,172800],
['vacuum or N₂ (thermal anneal), or solvent vapour','block-copolymer thin film microphase-separating','hotplate/oven or solvent-vapour chamber','chips to 300 mm wafers (directed self-assembly)'],
'PS-b-PMMA, PS-b-PDMS and high-χ copolymers|segregation needs χN > 10.5|neutral brush layer for perpendicular domains','thermal anneal 150–250 °C, minutes–hours|solvent-vapour anneal at room temperature','chain diffusion sets ordering time|film thickness commensurate with period',
'','period L₀ = 10–100 nm|guiding patterns give 2–4× density multiplication|defect annihilation kinetics',
'molar mass (period)|block volume fraction|annealing temperature/time|film thickness','period and feature size|morphology (spheres, cylinders, lamellae, gyroid)|orientation|defect density',[5,8,16,17,21,25],'solvent vapours',
'Pattern density from thermodynamics rather than exposure dose — very low energy per feature.',[5, 50],'ordering in minutes–hours');
e(113,'M 2 po co al pe L P I',[15,230],A,[600,86400],
['cleanroom air; vacuum filling common','liquid crystal between alignment-coated substrates','rubbing or photoalignment; capillary or vacuum filling','cm² cells to display-scale panels'],
'nematic mixtures (5CB, E7), cholesterics, reactive mesogens|polyimide or photoalignment layers|chiral dopant sets pitch','polyimide bake 180–230 °C|fill in the isotropic phase (clearing point 35–100 °C) and cool slowly','capillary filling time ∝ viscosity × length²/gap|flow-induced alignment defects',
'photoalignment with linearly polarised UV/blue light, 0.1–5 J/cm²|switching 1–10 V at ~1 kHz|pixel pitch down to ~1 µm for patterned alignment','cell gap 2–20 µm set by spacers|pretilt 1–5°|anchoring energy 10⁻⁵–10⁻³ J/m²',
'alignment pattern|cell gap|liquid-crystal mixture|voltage','director profile|retardation/phase|switching time|tunable resonance',[5,8,12,20,21],'low hazard; solvents for polyimide',
'Mature, scalable display technology reused for tunable metasurfaces; drive power is µW/cm².',[1000, 10000],'filling minutes; switching ms');
e(114,'N M 2 0 1 me sc po pe al lc L',[15,60],A,[60,14400],
['air at controlled humidity and dew point','receding meniscus of a colloidal suspension over a patterned template','evaporation and controlled contact-line motion','cm² templates; transfer printing afterwards'],
'particles 10 nm–10 µm at 0.01–1 wt%|surfactant tunes contact angle and stability','substrate held 5–30 K above dew point to control evaporation','evaporation-driven flow accumulates particles at the contact line|capillary force traps particles in template recesses',
'','receding contact angle 20–60°|contact-line speed 0.1–100 µm/s|trap geometry sets particle number and orientation',
'contact-line speed|contact angle|temperature|particle concentration','placement yield (> 95 % achievable)|position accuracy|particle orientation|array pitch',[5,8,9,12,21],'low hazard',
'Places single pre-synthesised nanoparticles deterministically with no vacuum or etching.',null,'mm²–cm² per hour');
e(115,'N M 1 2 3 po me sc ox pe al L P',[15,250],A,[60,172800],
['air, solvent vapour or vacuum','self-assembling species confined by topographic or chemical templates','self-assembly guided by pre-patterns','chips to wafers'],
'block copolymers, colloids, nanowires or molecules|chemical contrast or trench walls with tuned wetting|template from lithography, anodic alumina or imprint','as for the assembling species (ambient to 250 °C anneal)','commensurability between template and natural period|defects when mismatch exceeds ~10 %',
'','template period up to several times the natural period|registration to template within a few nm',
'template geometry|surface chemistry|commensurability|annealing','placement accuracy|long-range order|orientation|defect density',[5,8,16,17,21,25],'solvents and etchants from template making',
'Sparse templates multiply into dense patterns, cutting lithography effort; resolution is set jointly by template and assembling species.',[10, 100],'ordering minutes–hours');
e(116,'M N 3 1 co me ox po al pe lc L P',[15,80],A,[0.1,3600],
['air or sealed liquid cell','magnetically responsive colloids or composites in a carrier','static, rotating or gradient magnetic fields','mm cells to m-scale cured sheets'],
'magnetite or carbonyl-iron particles 10 nm–10 µm at 0.1–40 vol%|non-magnetic particles in ferrofluid (inverse systems)','ambient; curing of elastomer matrix 20–80 °C','chains, sheets and lattices form as field rises|rotating fields give planar assemblies|Mason number: viscous vs magnetic forces',
'flux density 1–500 mT|rotation 1–100 Hz|gradients for transport','cured magnetoactive elastomers change stiffness 10–100s of % under field|remanent programming of hard-magnetic composites',
'field strength and geometry|particle loading|rotation frequency|carrier viscosity','structure type|spacing|field-tunable stiffness or band gap|response time',[5,13,14,15,20,21],'strong magnetic fields',
'Remote, reversible structuring with abundant iron-based fillers.',null,'assembly ms–minutes');
e(117,'M N 2 3 po ox co pe al L',[15,40],A,[1,3600],
['sealed thin liquid cell','colloidal suspension between transparent electrodes','AC or DC electric field','cm² cells'],
'particles 0.1–10 µm at 0.1–30 vol%|low-ionic-strength water or organic solvent|conductivity 10⁻⁴–10⁻² S/m','ambient; Joule heating limits field in conductive media','electrohydrodynamic flows gather particles below ~1–10 kHz|dipolar chaining at higher frequency|electrophoretic compaction under DC',
'field 10³–10⁵ V/m|100 Hz–1 MHz|ITO electrode gap 10–200 µm','reversible: structure melts when the field is removed|electrode reactions above ~1 V DC',
'frequency|field amplitude|ionic strength|volume fraction','lattice type|lattice constant|switching time|optical stop band',[5,8,9,11,12,20],'low hazard',
'Reconfigurable photonic lattices at microwatt power; areas are limited.',null,'reconfiguration 0.1–100 s');
e(118,'N M 0 3 po co bi ox pe L P',[15,80],[1e5,5e5],[60,86400],
['closed microchannels','immiscible phases forming monodisperse droplets','syringe or pressure-driven pumps','single chips (mL/h) to parallelised L/h'],
'aqueous/oil phases with 0.5–5 wt% surfactant|monomer, polymer or nanoparticle payload|solidification by UV, gelation or solvent extraction','ambient to 80 °C for thermal curing','T-junction, flow-focusing or co-flow geometry|capillary number 10⁻³–10⁻¹ (dripping regime)|channels 10–500 µm',
'','channel wetting must favour the continuous phase|coefficient of variation < 1–3 %',
'flow-rate ratio|capillary number|channel size|surfactant','droplet/particle diameter 5–500 µm|monodispersity|core–shell structure|composition',[8,9,18,20,21,26],'fluorinated oils and surfactants',
'Exact dosing with almost no reagent waste; parallelisation is needed for meaningful output.',null,'droplets at 0.1–10 kHz');
e(119,'M N 3 ce co po ca al 3a lc L P I',[-196,0],[10,1.1e5],[3600,259200],
['air while freezing; vacuum < 600 Pa during sublimation','aqueous (or camphene, tert-butanol) slurry frozen directionally','cold finger or plate; freeze-dryer','cm samples to continuous tape freeze-casting'],
'solids 5–50 vol%|binder (PVA) and dispersant 0.5–2 wt%|additives (sucrose, ethanol, salts) change ice-crystal morphology','cold plate −5 to −196 °C; cooling 1–20 K/min|sublimation at −50 to −20 °C shelf, 10–100 Pa|sintering afterwards for ceramics','particles rejected by the ice front below a critical velocity, engulfed above|front velocity 1–100 µm/s',
'optional magnetic or electric alignment during freezing','lamellar wall spacing 5–100 µm, falling as front velocity rises|porosity 30–95 %|strong anisotropy',
'freezing-front velocity|solids loading|additives|temperature gradient','pore/lamella spacing|porosity|wall thickness|alignment',[5,9,10,14,16,21],'cryogens|vacuum',
'Water is the template and leaves no residue; freeze-drying energy is the main burden.',null,'front 1–100 µm/s');
e(120,'N M 3 po ca co ce bi 3a lc L P',[-196,-5],[10,1.1e5],[3600,259200],
['air while freezing; vacuum during drying','polymer solution, gel, emulsion or graphene-oxide dispersion with ice as template','freezers, cold baths, bidirectional cold fingers; freeze-dryer','vials to sheets and monoliths'],
'polymer or nanosheet 0.2–10 wt% (GO 2–10 mg/mL, cellulose nanofibres, chitosan)|cryo-polymerisation or crosslinking in the unfrozen phase|emulsion droplets add a second pore level','cryogelation −10 to −20 °C; rapid freezing in liquid N₂ gives fine pores|freeze-drying 10–100 Pa','solute concentrates between ice crystals|bidirectional gradients give long-range lamellae',
'','hierarchical pores from nm (walls) to 1–100 µm (ice)|aerogel density 2–100 kg/m³|elastic recovery of lamellar monoliths',
'freezing rate and direction|solute concentration|crosslinking|second template','pore-size hierarchy|density|compressive elasticity|conductivity',[5,8,9,10,14,21],'cryogens',
'Water-based route to ultralight foams; batch freeze-drying limits throughput.',null,'freeze minutes–hours; dry 12–72 h');
e(121,'N M 0 3 2 sc ox po bi pe hp L P I',[0,150],A,[600,6e5],
['air or inert; controlled vapour','solution, melt or colloidal fluid crystallising on seeds, templates or under gradients','controlled cooling, evaporation, antisolvent addition or fields','vials to industrial crystallisers'],
'supersaturation ratio 1.01–2 within the metastable zone|seeds 0.1–5 wt%|additives and templates select polymorph and face','cooling 0.05–1 K/min|temperature cycling dissolves fines','growth vs secondary nucleation balance|mixing avoids local supersaturation spikes',
'','template lattice match within a few % for epitaxy|attrition in stirred vessels',
'supersaturation|seeding|cooling/evaporation rate|additives','crystal size|polymorph|orientation|superlattice order',[3,6,7,16,21,26],'solvent hazards',
'Purifies as it structures; solvent recovery sets the footprint.',null,'growth 1 nm/s–1 µm/s');
e(122,'N M 3 ox po ca ce 3a L P I',[20,80],[1e5,1.5e7],[86400,1.2e6],
['air during gelation; supercritical CO₂ or ambient drying','wet gel network whose pore liquid is removed without collapse','sol–gel chemistry; high-pressure dryer','cm monoliths to m-scale blankets and granules'],
'silica (TEOS, waterglass), resorcinol–formaldehyde, cellulose, graphene gels|solvent exchange to ethanol or acetone, then CO₂|silylation (TMCS, HMDZ) for ambient-pressure drying','gelation and ageing 20–60 °C for 1–7 days|supercritical drying above 31.1 °C and 7.38 MPa (typically 40–50 °C, 10–15 MPa)','capillary pressure 2γ/r reaches tens of MPa in nm pores on evaporative drying|diffusion-limited solvent exchange ∝ thickness²|slow depressurisation',
'','porosity 90–99.8 %; density 3–500 kg/m³|pores 2–50 nm; surface area 500–1200 m²/g|brittle unless fibre-reinforced',
'precursor concentration|catalyst/pH|ageing time|drying route','density|pore size|surface area|thermal conductivity (≈ 0.012–0.02 W/m/K silica)',[5,7,8,9,10,21],'high-pressure CO₂|flammable solvents',
'Product saves energy in use (insulation), but solvent exchange is solvent-intensive; ambient-pressure drying and recycle loops cut impact.',null,'days per batch');
e(123,'N M 2 me po ox ca pe lc L P',[15,60],A,[10,14400],
['air','liquid–liquid or liquid–air interface','interfacial energy reduction; triggered by solvent or salt addition','cm² films; continuous transfer demonstrated'],
'nanoparticles 2–100 nm with tuned ligands|destabiliser (ethanol, salt) drives particles to the interface|interfacial polymerisation or 2-D polymer/COF formation','ambient','adsorption energy πr²γ(1 − |cos θ|)² ≫ kT for r > ~5 nm|Marangoni compaction of rafts',
'optional potential control at polarised liquid–liquid interfaces','interfacial tension 20–50 mN/m|film transfer by lifting or stamping|cracks on transfer',
'particle wettability (contact angle)|destabiliser dose|interfacial tension|particle concentration','packing density|interparticle gap|film area|layer number',[5,8,9,12,17,21],'organic solvents',
'Fast, equipment-free monolayers for plasmonic films and membranes.',null,'films form in seconds–minutes');
e(124,'N M 0 1 3 po bi co pe hp L P',[0,120],A,[1,604800],
['air or inert','solution of molecular building blocks','non-covalent recognition (host–guest, hydrogen bonding, π–π, metal–ligand)','mg to industrial (e.g. porous frameworks)'],
'building blocks µM–mM|association constants 10²–10¹² M⁻¹|stoichiometry and solvent polarity decide the assembly|modulators slow framework growth','ambient to 120 °C; annealing cycles correct errors','reversible bonding allows error correction|kinetic traps at high supersaturation',
'','mechanically soft; guest removal may collapse pores',
'building-block design|concentration|solvent|temperature','assembly size and topology|pore size|stimulus response|stability',[3,6,8,17,20,21],'solvent hazards (DMF)',
'Self-correcting and recyclable in principle; replace DMF with water or alcohols.',null,'seconds–days');

// ── H. High-pressure, vacuum, controlled-atmosphere and specialised ──
e(125,'N 0 1 3 ox ce sc F L P I',[100,450],[1e5,3.5e7],[3600,6e5],
['sealed liquid plus vapour; optional gas overpressure','pressurised aqueous or solvent solution','external heaters; magnetic-drive stirring','25 mL PTFE-lined bombs to multi-m³ vessels'],
'mineralisers (NaOH, KOH, NH₄F)|corrosion-compatible liner|gas-generating reactions raise pressure','PTFE liners ≤ 220–250 °C; alloy vessels to 400–500 °C|pressure follows the solvent isochore for the chosen fill','stirring 0–1000 rpm|temperature gradient drives dissolution–transport–growth (quartz)',
'','fill fraction 50–80 %; never liquid-full on heating|relief valve or burst disc|pressure-vessel inspection codes',
'temperature|fill fraction (pressure)|mineraliser|time','crystal size|phase|aspect ratio|yield',[6,7,9,10,21,27],'stored-energy vessel rupture|hot caustic release',
'Enables water-based crystal growth; heat integration and long vessel life offset capital and energy.',null,'batch 1–72 h');
e(126,'N 3 0 ca ce me hp L P I',[1200,1800],[4.5e9,8e9],[60,6e5],
['solid pressure medium; no free gas','solid or metal-solvent melt under quasi-hydrostatic pressure','belt, cubic or split-sphere press with internal graphite heater','mm³ research cells to cm³ industrial capsules'],
'graphite + Fe–Ni–Co solvent-catalyst for diamond; hBN + catalyst for cBN|nitrogen getters (Ti, Al) for colourless diamond|seed crystals','diamond 1300–1600 °C at 5–6 GPa|temperature gradient 20–50 K drives carbon transport to the seed|direct conversion needs ≥ 15 GPa and ~2300 °C','carbon diffuses through molten metal film|pressure medium (pyrophyllite, NaCl) transmits load',
'heater current kA at a few volts','tungsten-carbide anvils near their strength limit|pressure calibrated by fixed-point transitions|slow decompression avoids cracking',
'pressure|temperature|temperature gradient|catalyst composition','crystal size|nitrogen content/colour|grain size (polycrystalline)|inclusion content',[3,6,10,14,16,21],'anvil failure with large stored energy',
'Electricity for heating and pressing dominates; cobalt in catalysts and carbide tooling is a critical raw material.',null,'grit in minutes; gem growth ≈ mg/h over days');
e(127,'N M 2 0 sc me ox ep hp F L P',[-269,1200],[1e-9,1e-7],[3600,6e5],
['ultra-high vacuum; residual gas mainly H₂','atomically clean surfaces','turbomolecular, ion and sublimation pumps; bakeout heaters','research chambers to production MBE clusters'],
'contamination: one monolayer in ~1 s at 10⁻⁴ Pa, hours at 10⁻⁸ Pa|surface prepared by sputter–anneal cycles or flashing','bakeout 120–200 °C for 24–48 h|samples from 4 K (cryostats) to > 1200 °C (flash)','molecular flow; mean free path ≫ chamber size|outgassing and permeation set base pressure',
'in-situ electron and photon probes (LEED, XPS, RHEED, STM)','all-metal seals (CF flanges), low-outgassing materials|load-locks preserve vacuum',
'base pressure|bakeout|surface preparation|sample temperature','surface cleanliness|interface abruptness|impurity level|reconstruction',[3,7,17,21,23],'implosion and cryogen hazards|high voltage',
'A condition rather than a route; the temperature window spans cryostats to flash anneals. Continuous pumping and bakeouts are the energy cost.',null,'pump-down and bake 1–3 days');
e(128,'N 0 2 sc me po F hp L P',[15,40],[1e5,1.03e5],[60,6e5],
['Ar or N₂ with O₂ and H₂O < 1 ppm','dry inert gas at slight overpressure','circulation blower with catalyst and molecular-sieve purifier','single boxes to linked production lines'],
'handles alkali metals, organometallics, sulphides, perovskite and battery precursors|solvent vapours poison the purifier — use traps','box near ambient; internal hotplates and furnaces as fitted','gas recirculation tens of m³/h|antechamber pump–refill cycles (three or more)',
'','overpressure 100–1000 Pa|glove permeation is the main leak path',
'O₂ and H₂O level|purifier regeneration|antechamber cycling|solvent load','oxidation state|surface cleanliness|reproducibility|shelf life',[8,21,23,26,27],'asphyxiant gas|pyrophoric contents',
'Prefer N₂ where chemistry allows (argon is costlier and more energy-intensive); regenerate the purifier with dilute H₂.',null,'antechamber cycle 5–30 min');
e(129,'N M 2 3 po ox bi co F L P I',[10,90],A,[60,2.6e6],
['air at 5–98 % relative humidity','humid air over a film, powder or part','humidifier/dehumidifier with refrigeration and heaters','benchtop to walk-in chambers'],
'water activity drives hydrolysis, swelling and deliquescence|sol–gel, perovskite and hygroscopic-polymer processing','10–90 °C; humidity ±1–3 %; dew-point control','water sorption/desorption kinetics by diffusion|condensation when surfaces fall below dew point',
'','hygroscopic swelling strains|accelerated ageing (e.g. 85 °C / 85 % RH)',
'relative humidity|temperature|exposure time|air velocity','film morphology|water content|crystallinity|degradation rate',[6,8,9,21,22,23],'low hazard; microbial growth at high humidity',
'A condition applied to other routes; desiccant-wheel systems save energy over refrigeration at low dew points.',null,'stabilisation 10–60 min');
e(130,'N M 3 0 2 me po co bi F L P',[-269,-50],[1e-6,1.1e5],[60,86400],
['liquid or cold-gas N₂ or He; vacuum cryostats','cryogenic liquid, cold gas or cold solid contact','liquid N₂ (−196 °C), liquid He (−269 °C) or cryocoolers','dewars to industrial cryomills'],
'reactions and diffusion suppressed; oxidation limited under N₂|cryoprotectants for biological samples','cooling 1–10⁴ K/s (plunge freezing vitrifies water above ~10⁴ K/s)|thermal contraction 0.2–1.5 % to 77 K','film boiling limits heat transfer (Leidenfrost)|cryo-pumping of residual gases',
'','ductile-to-brittle transition aids milling of polymers and ductile metals|thermal-shock cracking',
'temperature|cooling rate|hold time|atmosphere','grain size (suppressed recovery)|amorphous/vitrified state|particle size|retained phases',[10,14,16,21,23],'cold burns|asphyxiation and pressure build-up',
'A condition applied to other routes; liquid-N₂ production costs ~0.5 kWh/kg, so recover cold and use closed-cycle coolers.',null,'cool-down seconds–hours');
e(131,'N 3 0 1 ce me ca co P I L',[400,2200],[10,1.5e5],[1800,172800],
['flowing Ar, N₂, NH₃, forming gas, CO/CO₂ or vacuum','solid powders or compacts in a gas-tight hot zone','tube or chamber furnace with graphite, Mo or W heaters','lab tubes to continuous pusher furnaces'],
'carbothermal reduction 1400–1600 °C|nitridation in N₂/NH₃; carburisation|oxygen and dew point at ppm level','ramp 1–20 K/min; hold 1–24 h|radiation dominates above ~1000 °C','gas flow 0.05–5 L/min; solid–gas reaction limited by pore diffusion|gaseous products (CO) swept away',
'','tube and seal limits: alumina ≤ 1700 °C; graphite higher in inert gas only|grain growth vs conversion',
'temperature|gas composition|flow rate|time','phase|particle/grain size|stoichiometry (C, N, O)|morphology (whiskers, powders)',[6,7,10,16,21,26],'toxic and flammable gases (CO, NH₃, H₂)|hot zone',
'High thermal demand; reactive precursors lower temperature and off-gas CO can be recovered as fuel.',null,'kg per batch to continuous t/day');
e(132,'N M 3 0 po ox co cf L P I',[32,150],[7.4e6,3.5e7],[600,86400],
['supercritical CO₂, optionally with 1–10 % co-solvent','supercritical fluid','high-pressure CO₂ pumps and heaters','mL cells to m³ extractors and dryers'],
'CO₂ is non-polar: dissolves fluorinated and small organic species; co-solvents extend range|rapid expansion (RESS) or antisolvent (SAS) precipitation|polymer plasticisation and foaming','above 31.1 °C and 7.38 MPa; typically 35–80 °C, 8–30 MPa|Joule–Thomson cooling on expansion','density 200–900 kg/m³ tunable with pressure|diffusivity ~10⁻⁸–10⁻⁷ m²/s; zero surface tension',
'','depressurisation rate controls foam cell size and gel cracking|CO₂ embrittlement of elastomer seals',
'pressure (density)|temperature|co-solvent|depressurisation rate','particle size|pore/cell size|residual solvent|surface area',[6,8,9,10,21,26],'high pressure|CO₂ asphyxiation',
'Replaces organic solvents with recyclable CO₂ (> 95 % recirculated in closed loops); compression energy is the cost.',null,'cycles 0.5–24 h');
e(133,'N M 2 3 me sc ox ce hp L P I',[150,1200],[1e-6,1e-1],[600,86400],
['high vacuum','solid film or part radiatively heated','resistance or lamp heaters in a pumped chamber','lab tubes to production vacuum furnaces'],
'avoids oxidation; desorbs hydrogen, water and organics|volatile elements (Zn, Mg, Mn) may evaporate|interdiffusion in multilayers','ramp 1–20 K/min; hold 0.5–24 h|radiative heat transfer only — slow, uniform cooling','outgassing load raises pressure during ramp|grain-boundary and interdiffusion govern microstructure',
'','stress relaxation and dewetting of thin films|thermal-mismatch stress on cooling',
'temperature|time|base pressure|ramp/cool rates','grain size|residual stress|interface sharpness|defect density',[10,14,16,21,22],'hot surfaces; implosion',
'Clean and gas-free; pump and heater electricity scale with cycle time, so batch loading matters.',null,'cycles 1–24 h');
e(134,'N 0 3 2 me ox ca sc hp L P I',[150,1100],[1e3,1e7],[600,86400],
['pure H₂ or forming gas (3–5 % H₂ in N₂/Ar)','reducing gas over powders, films or catalysts','furnace or pressurised reactor','lab tubes to industrial reduction furnaces'],
'reduces oxides of Cu, Ni, Co, Fe, W, Mo (Ellingham); not Al, Ti|dew point below −40 °C keeps H₂/H₂O ratio reducing|hydrogenation and defect engineering (e.g. black TiO₂)','150–1100 °C by oxide; low temperature preserves nanostructure|reduction of oxides is mildly endo- or exothermic','gas–solid reaction with water-vapour removal by flow|pore diffusion in powder beds',
'','hydrogen embrittlement of steels|leak-tight design, purge before and after',
'temperature|H₂ partial pressure/dew point|time|flow rate','degree of reduction|particle size|oxygen content|surface area',[6,7,9,10,21,27],'flammable 4–75 % in air; wide explosive range|high-pressure hydrogen',
'Water is the only by-product; climate benefit depends on low-carbon hydrogen.',null,'batch 0.5–24 h');
e(135,'N M 2 0 ox ce sc ep hp L P',[300,1500],[1e-4,1.1e5],[600,172800],
['O₂/inert mixtures, CO/CO₂ or H₂/H₂O buffers; vacuum with O₂ leak','oxide solid equilibrating with a gas of defined oxygen activity','furnace or deposition chamber with gas blending and zirconia oxygen sensor','lab reactors to production kilns'],
'oxygen partial pressure 10⁻²⁰–10⁵ Pa fixes oxidation state and vacancy concentration|defect chemistry: conductivity ∝ pO₂^(±1/n)|phase-stability windows (e.g. Fe₃O₄ vs Fe₂O₃)','300–1500 °C; equilibrium pO₂ of buffers varies strongly with temperature|quench preserves high-temperature state','oxygen surface exchange and chemical diffusion; equilibration time ∝ thickness²/D',
'zirconia sensor Nernst voltage E = (RT/4F) ln(pO₂/p_ref)','chemical expansion on reduction can crack ceramics',
'oxygen partial pressure|temperature|time|cooling path','oxygen stoichiometry|valence state|conductivity/magnetism|phase',[3,6,11,13,16,21],'toxic CO buffers|oxidiser-enriched atmospheres',
'Adds little energy beyond the furnace itself and prevents off-specification product.',null,'equilibration minutes–days');
e(136,'N M 2 3 0 me sc ce po ca L',[-196,800],[1e-6,1.1e5],[1,3.6e5],
['vacuum (ion and electron beams) or air/inert (gamma, X-ray)','solid under energetic particle or photon irradiation','ion accelerators, electron microscopes or accelerators, gamma sources, reactors','TEM foils to industrial e-beam and gamma processing'],
'displacement damage, implantation, radiolysis and cross-linking|radiolytic reduction of metal ions to nanoparticles in solution','−196 to 800 °C stage; temperature decides defect mobility|beam heating','radiation-enhanced diffusion|ballistic mixing at interfaces',
'ions 30 keV–MeV: 0.1–100 dpa at 10⁻⁴–10⁻² dpa/s|electrons 80–300 keV (TEM) to 10 MeV|gamma/e-beam 1–1000 kGy','swelling, hardening, amorphisation|in-situ TEM observation',
'particle type and energy|dose/fluence|dose rate|temperature','defect density|amorphous fraction|cross-link density|pore or bubble size',[16,17,21,22,23,27],'ionising radiation|activated samples',
'Reagent-free modification; shielding, licensing and source stewardship dominate cost and risk.',null,'seconds–days by dose');
e(137,'N 3 0 ca ce me sc ox L',[-269,5000],[1e8,4e11],[60,6e5],
['pressure medium He, Ne, Ar, or silicone oil/methanol–ethanol','microscopic sample between opposed diamond anvils','screw or gas-membrane loading; laser or resistive heating','sample chambers 10–300 µm; pL–nL volumes'],
'new high-pressure phases, hydrides and nanomaterials; many not recoverable to ambient|reaction with the medium or gasket possible','cryostat (4 K) to resistive heating (~1000 °C) to laser heating (> 5000 K); the window covers all modes|steep gradients in laser-heated spots','no flow; He and Ne stay quasi-hydrostatic to tens of GPa|chemical diffusion in laser-heated spots',
'optical, X-ray and Raman access through the anvils|ruby fluorescence or X-ray pressure standards','culets 30–500 µm (smaller gives higher pressure, routinely > 100 GPa, to ~400 GPa)|Re or steel gasket pre-indented to 20–40 µm|anvil failure near the limit',
'pressure|temperature|pressure medium|compression rate','crystal phase|density|grain size|recoverability',[3,6,14,16,17,23],'lasers and synchrotron radiation',
'Negligible material and energy use; a discovery tool rather than a production route.',null,'ng–µg per load');
e(138,'N M 3 0 po sc me bi co F L',[15,40],A,[5,1.6e7],
['cabin air (~101 kPa) or sealed experiment cells','fluids, colloids, melts and crystals under reduced gravity','orbital platforms, sounding rockets, parabolic flights, drop towers','mL-scale experiment cartridges'],
'protein, colloidal and semiconductor crystal growth|containerless (levitated) melts avoid crucible contamination','cabin 18–27 °C; furnaces and levitators reach > 1500 °C separately','buoyant convection and sedimentation suppressed: diffusion-dominated transport|Marangoni convection remains at free surfaces',
'optional electromagnetic or electrostatic levitation','residual acceleration ~10⁻⁶ g (orbit), 10⁻²–10⁻⁵ g (aircraft, drop tower)|g-jitter from crew and machinery',
'residual acceleration|experiment duration|temperature profile|concentration','crystal perfection and size|colloidal lattice order|compositional homogeneity',[8,9,16,18,23,24],'launch and containment requirements',
'Research environment with a very high launch footprint per gram; the time window spans drop towers (seconds) to orbit (months).',null,'seconds (drop tower) to months (orbit)');

// ── I. Natural and geological formation environments ──
e(139,'N 0 3 ox ce X',[700,1250],[1e5,1e9],[1,3e12],
['magmatic volatiles H₂O, CO₂, SO₂, HCl, HF','silicate melt, eruption plume and condensing vapour','geothermal heat; decompression','natural: m³ vents to km³ eruptions'],
'melt 45–77 wt% SiO₂|volatile exsolution on decompression|oxygen fugacity near the fayalite–magnetite–quartz buffer','rhyolite 700–900 °C to basalt 1100–1250 °C|plume quench up to 10³–10⁶ K/s forms glass and nanolites','melt viscosity 10–10¹² Pa·s|fragmentation when bubble fraction ≈ 75 % or strain rate exceeds relaxation|plume transport over 10–1000s km',
'volcanic lightning and plume charging','lithostatic pressure 0.1–1000 MPa|explosive decompression',
'magma composition|volatile content|decompression rate|cooling rate','ash particle size|glass vs crystal fraction|nanolite content|surface salts',[6,7,9,10,16,29],'respirable ash and crystalline silica|acid gases',
'Natural analogue for flame and plasma synthesis; not controllable. The time window spans plume quench (seconds) to magma storage (~10⁵ years).',null,'eruption ms–s; storage up to ~10⁵ yr');
e(140,'N 0 1 3 ox ce sc X',[50,600],[1e6,3e8],[3600,3e13],
['aqueous fluid with CO₂, H₂S, CH₄','hot brine in fractured and porous rock','geothermal gradient or magmatic heat','natural: veins to district-scale systems'],
'salinity 0–50 wt% NaCl-equivalent|pH 2–10|metals carried as chloride and bisulphide complexes|silica solubility falls on cooling','50–600 °C; boiling and cooling trigger precipitation|gradients 10–100 K/km regionally, far steeper near intrusions','Darcy flow in fractures; permeability 10⁻¹⁸–10⁻¹² m²|fluid mixing and depressurisation|water–rock reaction fronts',
'natural redox and self-potential gradients','hydrostatic to lithostatic pressure 1–300 MPa|hydraulic fracturing and vein sealing cycles',
'temperature|pressure|fluid chemistry|fluid–rock ratio','mineral phase|crystal size|vein texture|ore grade',[6,7,9,10,16,29],'toxic metals (As, Hg) and H₂S',
'Natural counterpart of hydrothermal reactors; geological timescales (to ~10⁶ years) are included in the window.',null,'colloid formation in seconds; veins over 10³–10⁶ yr');
e(141,'N 0 sc me ox X',[2,407],[2e7,4e7],[0.1,3e10],
['seawater mixing with vent fluid rich in H₂S, H₂, CH₄, CO₂','turbulent buoyant plume of hot acidic fluid in cold seawater','magmatic heat driving seawater circulation','natural: chimneys of metres; plumes over 1000s km'],
'vent fluid pH 2–4 (alkaline fields pH 9–11)|Fe, Mn, Cu, Zn at µmol–mmol/kg; H₂S 3–12 mmol/kg|reducing fluid meets oxygenated seawater','black smokers 350–400 °C entering 2 °C seawater|quench within seconds','turbulent entrainment dilutes ~10⁴-fold|pyrite, sphalerite and chalcopyrite nanoparticles resist settling and oxidise slowly',
'electrical potential across conductive chimney walls','hydrostatic pressure 20–40 MPa at 2000–4000 m depth|porous chimney walls',
'fluid temperature|metal and sulphide concentration|mixing rate|pH','particle size (many < 200 nm)|sulphide phase|aggregation|oceanic transport distance',[6,7,8,9,10,29],'not applicable (remote); toxic sulphide fluids',
'Natural flow reactor with extreme quench; informs continuous hydrothermal synthesis. Chimney lifetimes of years–centuries are included in the window.',null,'nucleation in < 1 s');
e(142,'N 2 0 1 ox ce X',[-10,40],[7e4,1e6],[3e7,3e13],
['air, soil gas with elevated CO₂, pore water','thin water films on mineral surfaces in soil and regolith','solar heating, rainfall, biological activity','natural: global land surface'],
'rain pH ~5.6; soil pH 3–9|organic acids and ligands from roots and microbes|incongruent dissolution: feldspar → kaolinite/smectite; Fe²⁺ → ferrihydrite (2–6 nm)','seasonal and diurnal cycling; freeze–thaw|rates roughly double per 10 K','unsaturated flow and wet–dry cycles|surface-reaction-limited dissolution 10⁻¹⁴–10⁻¹⁰ mol m⁻² s⁻¹|colloid mobilisation',
'UV photochemistry at exposed surfaces','frost, salt and root wedging create fresh surface',
'climate (temperature, rainfall)|parent mineralogy|pH and organic ligands|time','clay mineral type|particle size|surface area|iron-oxide phase',[6,7,8,9,22,29],'respirable mineral dust',
'Slow ambient synthesis of nanoclays; also the basis of enhanced-weathering carbon removal. Durations are years to millions of years.',null,'soil formation ~0.01–0.1 mm/yr');
e(143,'N 0 ox ca co X',[-80,40],[1e3,1.02e5],[600,1e6],
['air with trace SO₂, NH₃, amines, organics, iodine species','gas phase and suspended droplets','solar UV photochemistry; ionisation by cosmic rays','natural: boundary layer to stratosphere'],
'sulphuric acid 10⁶–10⁷ molecules/cm³ with NH₃, amines or oxidised organics|oxidants OH, O₃, NO₃|sea-spray and dust as primary particles','−80 to 40 °C; nucleation favoured when cold','new-particle formation 0.01–10 cm⁻³ s⁻¹ (higher in polluted or coastal air)|coagulation and condensation sinks|turbulent mixing and deposition',
'UV actinic flux|ion-induced nucleation','relative humidity controls droplet water content|residence time days–weeks',
'precursor vapour concentration|temperature|pre-existing aerosol surface|sunlight','particle number|size (1–100 nm on formation)|composition|cloud-condensation activity',[6,7,9,12,27,29],'fine-particle health burden',
'Natural aerosol reactor; guides exposure assessment and emission control for engineered nanomaterials.',null,'growth 1–20 nm/h');
e(144,'N 0 1 3 bi ox ce co X',[-2,40],[1e5,1e8],[600,3e9],
['body fluids, seawater or fresh water','confined biological compartments (vesicles, matrices, extracellular spaces)','metabolic energy (ion pumps, enzymes)','natural: nm crystals to reef-scale structures'],
'Ca²⁺, phosphate, carbonate, silicic acid and iron concentrated by transport|amorphous precursors (ACC, ACP)|acidic proteins and polysaccharides direct nucleation|pH ~7.4 (vertebrates), ~5 in silica deposition vesicles','organism temperature, −2 to 40 °C','membrane ion transport sets local supersaturation|matrix-confined diffusion',
'geomagnetic sensing by magnetite chains','organic matrix gives toughness far above the pure mineral|hydrostatic pressure to ~100 MPa in deep-sea organisms',
'ion supply|organic template|pH|confinement geometry','crystal size (bone apatite ≈ 50 × 25 × 3 nm; magnetosomes 35–120 nm)|orientation|polymorph|hierarchy',[3,6,7,8,14,29],'not applicable',
'Benchmark for ambient, water-based ceramic synthesis from abundant elements.',null,'diatom valve in hours; bone remodelling over years');
e(145,'N 0 ca co X',[300,1500],[7e4,1.02e5],[0.01,1e6],
['air with combustion products CO₂, CO, H₂O, organics','turbulent diffusion flame and smoulder zone','combustion of biomass','natural: wildfires of m² to 10⁴ km²'],
'fuel-rich pyrolysis forms polycyclic aromatics, then soot inception|black carbon, organic aerosol (tar balls), mineral ash (K, Ca salts)','flaming 800–1200 °C (peaks ~1500 °C); smouldering 300–700 °C','soot residence in flame ms–s|aggregation to fractal clusters (dimension ≈ 1.8)|plume lofting',
'thermal radiation dominated by soot','wind and fuel moisture control combustion efficiency',
'fuel type and moisture|oxygen supply|flame temperature|residence time','primary particle size 20–50 nm|aggregate size|organic/elemental carbon ratio|emission factor',[6,7,9,10,27,29],'smoke inhalation; carcinogenic polycyclic aromatics',
'Uncontrolled analogue of flame synthesis and a major climate and health forcing. The time window spans soot formation to multi-day fires.',null,'of order 0.5 g black carbon per kg dry fuel');
e(146,'N 0 2 3 ce ox X',[0,60],A,[86400,3e13],
['arid air; brine surface','evaporating natural brine','solar evaporation','natural: pans and sabkhas to basin-scale salt giants'],
'seawater sequence: carbonates → gypsum (≈ 3.5–4× concentration) → halite (≈ 10–11×) → K–Mg salts (> 60×)|ionic strength up to > 6 mol/kg; water activity ≈ 0.75 at halite saturation|microbial mats template mineral films','brine 0–60 °C with strong diurnal cycling','net evaporation ~1–3 m/yr in arid basins|density stratification and capillary wicking (efflorescence)',
'intense solar UV','crystallisation pressure in pores causes salt weathering|desiccation cracking',
'evaporation rate|brine composition|temperature|inflow/recharge','mineral sequence|crystal size and habit|nanoscale precursor phases|layer thickness',[6,7,8,9,16,29],'not applicable; saline dust',
'Solar-driven crystallisation — the model for solar-pond salt and lithium recovery. Durations span days (pans) to ~10⁶ years (basins).',null,'salt accumulation mm–cm per year');

// ── Condition categories (specification areas 30–34) ──
export const CONDITION_CATEGORIES = [
  { key: 'chemical', area: 30, title: 'Chemical and thermodynamic conditions', parameters: L('Chemical precursors|Concentrations and stoichiometric ratios|Chemical purity|Solvent selection|Solubility|Supersaturation|pH|Ionic strength|Chemical potential|Activity coefficients|Reaction equilibrium|Precipitation conditions|Dissolution conditions|Oxidation–reduction potential|Catalyst concentration|Surfactant and ligand concentrations|Dissolved gases|Reaction kinetics|Reaction enthalpy|Reaction duration|Feed composition and flow rate|Surface chemistry|Chemical contamination|Environmental reactivity') },
  { key: 'thermal', area: 31, title: 'Thermal and pressure conditions', parameters: L('Initial temperature|Operating temperature|Temperature gradients|Heating rates|Cooling rates|Maximum and minimum temperature|Thermal cycling|Heat flux|Heat-transfer coefficients|Joule heating|Laser heating|Microwave heating|Reaction-generated heat|Thermal radiation|Ambient temperature|Pressure|Pressure gradients|Hydrostatic pressure|Vapour pressure|Gas partial pressures|Vacuum levels|Compressibility|High-pressure operation|Pressure transients|Temperature and pressure safety limits') },
  { key: 'transport', area: 32, title: 'Fluid, porous-media and mass-transfer conditions', parameters: L('Fluid density|Viscosity|Surface tension|Volumetric flow rate|Fluid velocity|Reynolds number|Flow regime|Mixing intensity|Residence time|Fluid recirculation|Diffusion coefficients|Mass-transfer coefficients|Reaction–transport coupling|Permeability|Porosity|Tortuosity|Capillary pressure|Pore-fluid saturation|Concentration gradients|Precursor depletion|Particle sedimentation|Particle suspension|Aggregation|Gas evolution|Multiphase interfaces|Pressure drop') },
  { key: 'fields', area: 33, title: 'Electrical, electromagnetic and optical conditions', parameters: L('Electric-field magnitude and direction|Voltage|Current|Electrical frequency|Waveform|Electrode geometry|Electrical conductivity|Electric potential|Ionic mobility|Charge distribution|Magnetic-field magnitude and direction|Magnetic bias|Electromagnetic wavelength|Incident wave angle|Polarisation|Optical intensity|Illumination duration|Pulsed versus continuous excitation|Resonance frequency|Electromagnetic power|Dielectric heating|Electromagnetic shielding|Optical damage thresholds|Frequency-dependent losses') },
  { key: 'mechanical', area: 34, title: 'Mechanical, atmospheric, biological and service conditions', parameters: L('Tensile and compressive loads|Shear and torsion|Bending|Pre-strain|Dynamic loading|Vibration frequency|Impact loading|Mechanical constraints|Contact pressure|Friction|Fatigue cycles|Ambient humidity|Oxygen concentration|Inert-gas environment|Corrosive atmosphere|Radiation exposure|UV exposure|Biological fluids|Cell-culture conditions|Microbial exposure|Cleanroom contamination|Substrate condition|Storage environment|Operational duty cycle|Service duration|Long-term ageing') }
];

// ── Table J: metamaterial class → concrete creation environments ──
export const METAMATERIAL_ROUTES = [
  { cls: 'Electromagnetic metamaterials', envIds: [62, 61, 71, 47, 91, 90, 57, 33, 32, 93],
    note: 'Lithography (62, in cleanrooms 61); PCB-style fabrication = photolithography + wet copper etching (71) + electroplating (47); precision printing (91, 90, 57); thin-film deposition (33, 32); roll-to-roll (93) for large areas.' },
  { cls: 'Optical metamaterials', envIds: [63, 64, 32, 33, 29, 69, 54, 60],
    note: 'Electron-beam (63) and nanoimprint (64) lithography; deposition by evaporation (32), sputtering (33) and ALD (29); pattern transfer by RIE/ICP etching (69, 54); FIB (60) for prototypes.' },
  { cls: 'Photonic metamaterials', envIds: [61, 62, 65, 34, 27, 109, 112, 66, 67],
    note: 'Nanofabrication cleanrooms (61, 62, 65); epitaxial deposition by MBE (34) and CVD/MOCVD (27); self-assembly of colloidal crystals (109) and block copolymers (112); interference lithography (66) and two-photon writing (67) for 3-D photonic crystals.' },
  { cls: 'Acoustic metamaterials', envIds: [86, 87, 88, 85, 83, 96, 70, 94],
    note: 'Additive manufacturing in polymer (86, 87, 88, 85) and metal (83); micromachining (96, 70 for MEMS-scale devices); precision assembly of resonators and membranes (94).' },
  { cls: 'Mechanical metamaterials', envIds: [67, 68, 87, 86, 85, 83, 84, 89, 70, 96],
    note: '3-D printing (87, 86, 85, 83, 84, 89); microfabrication (70); laser writing (67 for micro/nanolattices, 68); precision manufacturing (96).' },
  { cls: 'Auxetic metamaterials', envIds: [88, 86, 87, 85, 83, 84, 95, 96],
    note: 'Polymer printing (88, 86, 87, 85); metal additive manufacturing (83, 84); lattice and sheet fabrication by kirigami cutting (95) and machining (96).' },
  { cls: 'Thermal metamaterials', envIds: [89, 92, 33, 29, 30, 44, 119, 94, 96],
    note: 'Multimaterial fabrication (89, 92, 94, 96); thin-film deposition (33, 29, 30); composite manufacturing by sintering (44) and freeze-casting (119).' },
  { cls: 'Seismic metamaterials', envIds: [88, 96, 94],
    note: 'Field-scale geotechnical construction, drilling, foundation and periodic-barrier installation lie outside the 146 listed environments; the ids given are the routes for laboratory scale models (printed or machined resonator arrays and assembled inclusions).' },
  { cls: 'Elastic-wave metamaterials', envIds: [96, 83, 86, 87, 88, 94, 70],
    note: 'Precision machining (96); additive manufacturing (83, 86, 87, 88); embedded resonator assembly (94); deep reactive ion etching (70) for chip-scale phononic crystals.' },
  { cls: 'Programmable metamaterials', envIds: [94, 89, 87, 91, 95, 116, 113],
    note: 'Robotics and microassembly (94); multimaterial printing (89, 87, 91); smart-material integration through origami/kirigami with active hinges (95), magnetically programmed composites (116) and liquid-crystal cells (113).' },
  { cls: 'Active and tunable metamaterials', envIds: [61, 62, 75, 76, 113, 116, 87, 91],
    note: 'Semiconductor fabrication (61, 62, 75, 76) for embedded diodes, MEMS and transistors; printed embedded electronics (91); responsive polymer and liquid-crystal processing (113, 116, 87).' },
  { cls: 'Nonlinear metamaterials', envIds: [63, 62, 69, 34, 35, 76, 94, 95],
    note: 'Micro/nanofabrication (63, 62, 69); nonlinear material integration by epitaxy (34), pulsed laser deposition (35) and wafer bonding (76); architected mechanical assembly of bistable elements (94, 95).' },
  { cls: 'Chiral metamaterials', envIds: [67, 68, 63, 72, 94, 103],
    note: '3-D direct laser writing (67, 68); lithography (63) with multilayer alignment; focused electron-beam deposition of nanohelices (72); precision assembly (94) and DNA-origami placement of plasmonic particles (103).' },
  { cls: 'Topological mechanical and acoustic metamaterials', envIds: [87, 86, 88, 83, 96, 70, 67],
    note: 'Patterned lattice manufacturing by machining (96) and deep etching (70); additive fabrication (87, 86, 88, 83) and two-photon writing (67) at the microscale.' },
  { cls: 'Metasurfaces', envIds: [63, 62, 64, 65, 29, 32, 33, 69, 54, 91, 90, 93],
    note: 'Lithography (63, 62, 65); deposition (29, 32, 33); etching (69, 54); nanoimprinting (64, roll-to-roll 93); printed electronics (91, 90) for microwave and terahertz designs.' }
];

// ── Rule-based, explainable recommender ──
const KIND = { nanomaterial: 'N', metamaterial: 'M' };
const FAM = { metal: 'me', oxide: 'ox', semiconductor: 'sc', carbon: 'ca', polymer: 'po', ceramic: 'ce', composite: 'co', biological: 'bi' };
const NEED = { alignment: 'al', conformal: 'cf', periodic: 'pe', 'low-cost': 'lc', '3d-architecture': '3a', 'high-purity': 'hp', epitaxial: 'ep' };
const SCALE = { lab: 'L', pilot: 'P', industrial: 'I' };
const DIMS = ['0', '1', '2', '3'];
const DIM_WORDS = {
  '0D': ['nanoparticle', 'quantum dot', 'nanocrystal', 'nanopowder', 'particles'],
  '1D': ['nanorod', 'nanowire', 'nanotube', 'nanofibre', 'fibre'],
  '2D': ['film', 'graphene', 'monolayer', 'membrane', 'coating', 'metasurface', 'multilayer', 'nanolayer'],
  '3D': ['lattice', '3d', 'porous', 'aerogel', 'cellular', 'architect', 'scaffold', 'dense', 'bulk']
};
const FAM_WORDS = {
  metal: ['metal', 'gold', 'silver', 'copper', 'alloy'], oxide: ['oxide', 'zno', 'tio', 'silica', 'titania'],
  semiconductor: ['semiconductor', 'quantum'], carbon: ['carbon', 'graphene', 'diamond', 'fullerene', 'soot'],
  polymer: ['polymer'], ceramic: ['ceramic'], composite: ['composite'], biological: ['bio', 'protein', 'dna', 'peptide', 'lipid', 'tissue']
};
const KIND_WORDS = { metamaterial: ['metamaterial', 'metasurface', 'lattice', 'photonic', 'auxetic'], nanomaterial: ['nano', 'quantum'] };
const fmtNm = v => (v >= 1e6 ? `${+(v / 1e6).toPrecision(3)} mm` : v >= 1e3 ? `${+(v / 1e3).toPrecision(3)} µm` : `${+v.toPrecision(3)} nm`);
const VACUUM_FREE_PA = 5e4;

// Rank environments for a target; additive scoring, every contribution listed in why. All query fields optional.
// envList (spec.json environments) adds name/makes keyword evidence. Returns the top 12 [{id, score, why, name?}].
export function recommendEnvironments(q = {}, envList) {
  const { kind, dimensionality, materialFamily, scale } = q || {};
  const feat = Number(q?.featureSize_nm), Tmax = Number(q?.T_max_C);
  const needs = Array.isArray(q?.needs) ? q.needs.map(n => String(n).toLowerCase()) : [];
  const byId = {};
  if (Array.isArray(envList)) for (const x of envList) if (x && x.id != null) byId[x.id] = x;
  const out = [];
  for (const key of Object.keys(ENV_CONDITIONS)) {
    const id = Number(key), c = ENV_CONDITIONS[id], tags = ENV_TAGS[id].split(' ');
    const has = t => tags.includes(t);
    const dims = tags.filter(t => DIMS.includes(t)), fams = tags.filter(t => Object.values(FAM).includes(t));
    const info = byId[id], text = info ? `${info.name || ''} ${info.makes || ''}`.toLowerCase() : '';
    let score = 0; const why = [];
    const add = (pts, msg) => { score += pts; why.push(`${pts >= 0 ? '+' : '−'}${Math.abs(pts)} ${msg}`); };

    if (has('X')) add(-8, 'natural formation environment — a reference analogue, not an engineered production route');
    if (has('F')) add(-1, 'enabling facility or condition, normally combined with a primary synthesis route');

    if (KIND[kind]) {
      if (has(KIND[kind])) add(3, `established route for ${kind}s`);
      else add(-3, `not a usual route for ${kind}s`);
      if (text && KIND_WORDS[kind].some(w => text.includes(w))) add(1, `specification lists ${kind}-type products for this environment`);
    }
    if (DIM_WORDS[dimensionality]) {
      const d = dimensionality[0];
      if (dims.includes(d)) {
        add(3, `produces ${dimensionality} structures`);
        if (dims[0] === d) add(1, `${dimensionality} structures are its principal product form`);
      } else add(-2, `${dimensionality} structures are not typical products`);
      if (text && DIM_WORDS[dimensionality].some(w => text.includes(w))) add(1, `specification names ${dimensionality} products (${info.makes})`);
    }
    if (FAM[materialFamily]) {
      const f = FAM[materialFamily];
      if (fams.includes(f)) {
        add(3, `compatible with ${materialFamily} materials`);
        if (fams[0] === f) add(1, `${materialFamily} is its principal material family`);
      } else add(-2, `${materialFamily} materials are not typical for this environment`);
      if (text && FAM_WORDS[materialFamily].some(w => text.includes(w))) add(1, `specification names ${materialFamily} products`);
    }
    if (Number.isFinite(feat) && feat > 0 && c.resolution_nm) {
      const [lo, hi] = c.resolution_nm;
      if (feat < lo) add(-5, `required ${fmtNm(feat)} features are below its best resolution (${fmtNm(lo)})`);
      else if (feat < hi) add(2, `${fmtNm(feat)} features are achievable but near the resolution limit (${fmtNm(lo)}–${fmtNm(hi)})`);
      else if (feat > 1000 * hi) add(-1, `resolution (${fmtNm(lo)}–${fmtNm(hi)}) is far finer than the ${fmtNm(feat)} needed — slower and costlier than necessary`);
      else add(3, `resolves ${fmtNm(feat)} features comfortably (minimum ${fmtNm(lo)}–${fmtNm(hi)})`);
    }
    if (Number.isFinite(Tmax) && q.T_max_C !== null && q.T_max_C !== undefined && q.T_max_C !== '') {
      if (c.T_C[0] > Tmax) add(-5, `minimum process temperature ${c.T_C[0]} °C exceeds the ${Tmax} °C thermal budget`);
      else if (c.T_C[1] <= Tmax) add(2, `whole temperature window (${c.T_C[0]}–${c.T_C[1]} °C) fits the ${Tmax} °C thermal budget`);
      else add(1, `can run within the ${Tmax} °C thermal budget (lower part of its ${c.T_C[0]}–${c.T_C[1]} °C window)`);
    }
    if (SCALE[scale]) {
      if (has(SCALE[scale])) add(2, `practised at ${scale} scale (${c.scale})`);
      else add(-2, `not normally run at ${scale} scale (${c.scale})`);
    }
    for (const n of needs) {
      if (n === 'vacuum-free') {
        if (c.P_Pa[0] >= VACUUM_FREE_PA) add(2, 'runs at or above ambient pressure — no vacuum system needed');
        else add(-3, `needs vacuum (down to ${c.P_Pa[0].toExponential(0)} Pa)`);
      } else if (NEED[n]) {
        if (has(NEED[n])) add(2, `meets the need: ${n}`);
      }
    }
    const r = { id, score, why };
    if (info && info.name) r.name = info.name;
    out.push(r);
  }
  out.sort((a, b) => b.score - a.score || a.id - b.id);
  return out.slice(0, 12);
}
