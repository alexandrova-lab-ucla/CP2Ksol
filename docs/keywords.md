# Keyword reference

[English](keywords.md) · [中文](keywords.zh.md) · [Back to the README](../README.md)

This page lists every input keyword that the implementation adds or changes. Defaults are shown in the units the input uses. Aliases are in parentheses. Defaults and ranges are those of `src/input_cp2k_dft.F`, `src/input_cp2k_scf.F`, `src/input_cp2k_print_dft.F`, `src/start/input_cp2k.F` and the checks in `src/cp_control_utils.F`, `src/scf_control_types.F`, `src/qs_sccs.F`. The `&SCCS &DEVEL` section holds developer-only keywords; leave it at its defaults.

### `&DFT &SCCS`

#### Model and solver

| Keyword | Default | Meaning |
|---|---|---|
| `METHOD` | `ANDREUSSI` | Smoothing of the dielectric function. Two values added: `ISODENSITY`, `SOLVENT_EXCLUDED` |
| `SOLVER` | `JACOBI`; `NEWTON_CG` for `ISODENSITY` and `SOLVENT_EXCLUDED` when not given | `JACOBI`: damped fixed point on the polarisation charge, convergence test on its change per step (`EPS_SCCS`). `NEWTON_CG`: Newton on the true residual, inner CG, line search (`SOLVER_TOLERANCE`, `SOLVER_MAX_ITER`, `SCF_TOL_FACTOR`). `NEWTON_CG` needs `RELATIVE_PERMITTIVITY` > 1 and `METHOD` other than `SAA_ANDREUSSI` |
| `SOLVER_GRID_LEVEL` | 1 | MGRID level on which the solver works; 1 is the auxiliary basis grid. n > 1 restricts the solute density to level n and prolongates the solvation potential back; the vacuum Hartree energy stays on the fine grid. `NEWTON_CG` only. `SOLVENT_EXCLUDED` aborts for n > 1, `ISODENSITY` warns that convergence is poor |
| `SOLVER_STRIDE` | 1 | Solve on every n-th activated SCF step; in between the potential is frozen. Must be ≥ 1 |
| `PSI_FORM` | F; when not given T for `ISODENSITY` and for `SOLVENT_EXCLUDED` with `NONLINEAR_DIELECTRIC F`, with `NEWTON_CG` and without `&EXCLUSION` | ψ = √ε φ discretisation (https://doi.org/10.1063/1.4865107; https://doi.org/10.1063/1.4939125): the operator s(−∇²+q)s with q = ∇²(√ε)/√ε defines the discrete −∇·(ε∇), so the inner CG costs 2 FFTs per iteration. `NEWTON_CG` only; `NONLINEAR_DIELECTRIC` must be off. `PSI_FORM F` selects the plain operator |
| `SOLVER_TOLERANCE` | 1.0E-6 | RMS residual tolerance [e/bohr³] of the solver on a cold start and when `SCF_TOL_FACTOR` is 0 |
| `SOLVER_MAX_ITER` | 50 | Newton iteration cap |
| `SCF_TOL_FACTOR` | 0.1; when not given 0.3 for a linear dielectric without nonlinear ions | From the second activated SCF step the solver stops at RMS ≤ max(`SOLVER_TOLERANCE`, factor × RMS of the entry residual). 0 switches the coupling off. An explicit value always wins |
| `SOLVER_ENERGY_TOL` | 1.0E-10 Ha | Linear-dielectric / linearised-ion path: the CG also stops when its last 20 iterations raised the concave free-energy functional by less than this in total. 0 = RMS test only. Must be ≥ 0 |
| `SOLVER_MAX_FAILURES` | 0 | Abort after this many consecutive SCF steps whose solve did not converge; 0 = warn only. A step with a frozen potential neither counts nor resets |
| `GAMMA` (`SURFACE_TENSION`) | 0 mN/m | CP2K keyword; see [Surface tension](#surface-tension) |

#### Keywords inherited from CP2K

| Keyword | Default | Meaning |
|---|---|---|
| `RELATIVE_PERMITTIVITY` (`DIELECTRIC_CONSTANT`, `EPSILON_RELATIVE`, `EPSILON_SOLVENT`) | 80 | Bulk relative permittivity ε_b of the solvent. Use 78.4 for water at 298 K |
| `EPS_SCCS` (`EPS_ITER`, `TAU_POL`) | 1.0E-6 | Convergence tolerance of the `JACOBI` iteration on the polarisation density |
| `EPS_SCF` | 0.5 | The SCCS cycle starts once the SCF has converged to this threshold |
| `MAX_ITER` | 100 | Iteration cap of the `JACOBI` solver |
| `MIXING` (`ETA`) | 0.6 | Damping of the `JACOBI` iteration |
| `ALPHA` | 0 mN/m | Repulsion term αS |
| `BETA` | 0 GPa | Dispersion term βV |
| `DELTA_RHO` | 2.0E-5 | Numerical increment for the quantum surface of the cavity |
| `DERIVATIVE_METHOD` | `FFT` | Numerical derivatives on the real-space grid: `FFT`, `CD3`, `CD5`, `CD7` |

`METHOD` values `ANDREUSSI`, `FATTEBERT-GYGI` and `SAA_ANDREUSSI` and their sections are those of CP2K; see the CP2K manual.

#### Electrolyte

| Keyword | Default | Meaning |
|---|---|---|
| `DEBYE_LENGTH` (`LAMBDA_D`, `DEBYE_SCREENING_LENGTH`) | 0 Å = off | Linearised Poisson–Boltzmann screening length. 1 M 1:1 aqueous at 298 K: 3.04 Å. Needs `POISSON_SOLVER PERIODIC`, `RELATIVE_PERMITTIVITY` > 1, and `METHOD` other than `SAA_ANDREUSSI` |
| `ION_CONCENTRATION` | 0 = off | Bulk concentration in mol/l; when > 0 it replaces `DEBYE_LENGTH`, which is then λ_D = √(ε_b k_B T/(8π c z²)) |
| `ION_VALENCE` | 1 | Valence z of the symmetric z:z electrolyte (`NONLINEAR_IONS` and `ION_CONCENTRATION`) |
| `TEMPERATURE` | 298 K | Electrolyte temperature in the Boltzmann factors (`NONLINEAR_IONS`) and in λ_D from `ION_CONCENTRATION` |
| `NONLINEAR_IONS` (`NLPB`, `SIZE_MODIFIED_PB`) | F | Size-modified nonlinear ions instead of the linearisation. Needs `DEBYE_LENGTH` > 0 |
| `ION_STERIC_RADIUS` (`ION_RADIUS`) | 0 Å | Hydrated-ion radius of the lattice gas; 0 gives point ions (plain sinh nonlinear PB) |
| `ION_PACKING_DIAMETER` | −1 = 2^(5/6) × `ION_STERIC_RADIUS` | Close-packed ion diameter d_ion; n_max = 1/d_ion³. The input aborts when θ_b ≥ 1 |

#### Nonlinear dielectric

| Keyword | Default | Meaning |
|---|---|---|
| `NONLINEAR_DIELECTRIC` (`NLDIEL`, `DIELECTRIC_SATURATION`) | F; T for `SOLVENT_EXCLUDED` unless given | Rotational saturation of the dielectric. Needs `SOLVER NEWTON_CG` or `NONLINEAR_IONS`. Not combinable with `PSI_FORM` |
| `EPS_INF` | 1.78 | High-frequency relative permittivity |
| `SOLVENT_DENSITY` | 0.0335 Å⁻³ (4.9649E-3 bohr⁻³) | Solvent molecule number density |
| `SOLVENT_DIPOLE` | 0.9449 e·bohr | Solvent molecular dipole |
| `SMOOTHING_WIDTH` | 0.125 Å | Gaussian smoothing of the polarisation and bound charge; 0 = off |

#### Pseudo-core cavity density

| Keyword | Default | Meaning |
|---|---|---|
| `ADD_CORE_DENSITY` (`PSEUDO_CORE`) | T | Cavity density = valence + pseudo-core density. `ISODENSITY` aborts with `F` or `CORE_DENSITY_SCALE` ≤ 0 |
| `CORE_DENSITY_SCALE` (`PSEUDO_CORE_SCALE`) | 1 | Scales Z − Z_eff |
| `CORE_DENSITY_ZMIN` (`PSEUDO_CORE_ZMIN`) | 1 | Kinds with a smaller atomic number get no core |
| `CORE_DENSITY_WIDTH` (`PSEUDO_CORE_SIGMA`) | 0 = GTH core width | Gaussian width σ of the core density |
| `CORE_DENSITY_WIDTH_SCALE` (`PSEUDO_CORE_SIGMA_SCALE`) | 1 | Multiplies the GTH core width |
| `CORE_DENSITY_CHARGE` (`PSEUDO_CORE_CHARGE`) | −1 = `CORE_DENSITY_SCALE` × (Z − Z_eff) | Electrons per atom in the core Gaussian |

#### Speed-ups of the Poisson–Boltzmann solve

All are on by default except `NEWTON_CG_FORCING`; switching them off selects the plain code path. They change the CG path, so energies move at the level of the SCF stopping criterion.

| Keyword | Default | Meaning |
|---|---|---|
| `COARSE_DIRECT_SUM` | T | Two-level preconditioner of the ψ-form CG (coarse space: low plane waves s·exp(iGr)): projection and expansion as separable phase sums on the local grid block instead of 3-D FFTs, 4 instead of 6 FFTs per CG iteration. The phase convention is checked in the first solve and the FFT route is used if the check fails |
| `LAP_CARRY` | T | ψ-form: ∇²(sφ) after a CG run and at the trial points of the line search is updated from terms the CG already computed, which saves 2 FFTs per residual |
| `LINEAR_SINGLE_CG` | T | ψ-form, linear dielectric, no ions: one CG straight to the SCF-coupled tolerance. `F` reduces the residual by a factor of 10 per Newton step |
| `HISTORY_GUESS` | 2 | Linear dielectric, no ions: number (0–3) of previous potential increments of the same SCF run that span the initial guess; used only if it lowers the entry residual below 0.9 of its value. 0 = off |
| `LINEAR_CG_FLOOR` | 0.8 | Linear dielectric, no ions: the CG stops at RMS ≤ max(0.1 × entry RMS, value × tolerance). Range 0–1; 0 = stop at 0.1 × entry RMS |
| `COARSE_SPACE_IONS` | T | ψ-form with linearised Debye ions: use the two-level preconditioner; a CG iteration costs 4 FFTs instead of 3 and the iteration count drops by about half when the ion-accessible region fills a small part of the cell. Acts only when N_e can change during the SCF (`&SCF &GCE` with `MAX_OUTER_UPDATES` > 0, `UPDATE_EVERY SCF_STEP`, scheme 1 or `EARLY_UPDATE_DELTA` > 0), the same condition as `SCF_TOL_JUMP`; a fixed-N_e or non-GCE run always uses the plain preconditioner |
| `SCF_TOL_JUMP` | 100 | The SCF-coupled tolerance follows the minimum of the entry residuals since the last reset. An entry residual larger than this factor times the minimum restarts the minimum. 0 = off, otherwise > 1. Acts only with Debye ions and when `&SCF &GCE` can change N_e during the SCF (`MAX_OUTER_UPDATES` > 0, `UPDATE_EVERY SCF_STEP`, scheme `MU_MIXING`, or `EARLY_UPDATE_DELTA` > 0); an ion-free or fixed-N_e run behaves as 0 |
| `COLD_START_FROM_LAST` | T | The first solve of an SCF run starts from the potential at the end of the previous run; its tolerance stays `SOLVER_TOLERANCE` |
| `NEWTON_CG_FORCING` | 0 = off | Nonlinear paths: the CG of each Newton step stops at RMS ≤ max(0.1 × entry RMS, value × tolerance). Range 0–1. About 15 % faster, energies move by up to about 1E-5 Ha for `EPS_SCF` 1E-6 |

#### Diagnostics

| Keyword | Default | Meaning |
|---|---|---|
| `PRINT_POL_DIAGNOSTICS` | F | Print the terms of the polarisation energy and one line per CG solve; the energy terms are also printed whenever the positive-polarisation-energy warning fires |
| `SOLVER_DEBUG` | F | Finite-difference check of the Jacobian and free energy on the first Newton step of every solve |

### `&SCCS &ISODENSITY`

| Keyword | Default | Meaning |
|---|---|---|
| `RHO_CAVITY` (`NC`, `RHO_CUTOFF`) | 0.0030 Å⁻³ (4.44554E-4 bohr⁻³) | Density isovalue of the cavity boundary; fitted on GPW/GTH densities to 11 experimental solvation free energies |
| `SWITCH_WIDTH` (`SIGMA`) | 0.6 | Width of the erfc switch on the ln(n/`RHO_CAVITY`) axis |

### `&SCCS &SOLVENT_EXCLUDED`

| Keyword | Default | Meaning |
|---|---|---|
| `RHO_CAVITY` (`NC`) | 0.0158 Å⁻³ (2.34132E-3 bohr⁻³) | Density isovalue of the innermost surface |
| `SWITCH_WIDTH` (`SIGMA`) | 0.6 | Width of the erfc switch, used for every surface |
| `PROBE_KERNEL_WIDTH` (`KERNEL_A`) | 0.125 Å | Decay length a of the probe-shell kernel; the shell thickness is a/`SWITCH_WIDTH` |
| `SOLVENT_RADIUS` | 1.4 Å | Closest approach of a solvent molecule centre to the van der Waals surface |
| `ION_RADIUS` | 4.0 Å | Closest approach of an electrolyte ion; between this surface and the dielectric surface there is solvent but no ion |
| `DIELECTRIC_RADIUS` | 1.0 Å (0.7 Å in SE) | Inward displacement of the surface where dielectric screening starts; sets the width of the unscreened gap at an electrode. A linear dielectric with a value > 0 warns that the SCF may collapse |
| `CAVITATION_RADIUS` | 0 Å | Inward displacement of the surface that carries the cavitation term |

### `&SCCS &REGION`

Planar confinement of the implicit solvent for cells that contain explicit water. The dielectric and ionic cavities are multiplied by m(z) = ½[1 + erf((z − LOWER_EDGE)/(WIDTH√2))] × ½[1 + erf((UPPER_EDGE − z)/(WIDTH√2))] along the cell axis `AXIS`, summed over periodic images. `LOWER_EDGE`, `UPPER_EDGE` and `WIDTH` are lengths along the lattice vector of that axis, so the vector must be perpendicular to the other two. The section is active when it is present. Every SCCS call reports the cavity volume that lies outside the region and, after the mask is applied, readings in the window, in the whole cell and in the spheres around `PROBE_ATOMS`. The cavitation surface is not masked. Not available with `SAA_ANDREUSSI`.

| Keyword | Default | Meaning |
|---|---|---|
| `AXIS` | 3 | Cell axis (1, 2 or 3) normal to the explicit/implicit boundary |
| `LOWER_EDGE` | none | Lower edge of the region, in [0, L). Required unless `&FLEXIBLE` is present; give both edges or neither |
| `UPPER_EDGE` | none | Upper edge, > `LOWER_EDGE`; may exceed L so that the region wraps through the boundary; `UPPER_EDGE` − `LOWER_EDGE` < L − 6 `WIDTH` |
| `WIDTH` | 1.0 Å | Width of the erf edges; must be > 0. The film uses 0.5 Å |
| `APPLY` | T | `F` leaves the model unchanged and only reports the leakage |
| `PROBE_ATOMS` | none | Atoms (1-based) around which the readings are also taken; reporting only |
| `PROBE_RADIUS` | 3.0 Å | Radius of those spheres; must be > 0 |

#### `&REGION &FLEXIBLE`

Region boundary that follows the explicit solvent during a geometry optimisation or MD. The listed atoms define a coarse-grained density ρ̄(r) = Σᵢ wᵢ g(r − Rᵢ), with g a normalised Gaussian of width `COARSE_GRAINING_LENGTH`, summed over all periodic images. With t = erfc⁻¹(2ρ̄/`DENSITY`) the mask is m(r) = ½ erfc[(`OFFSET` − √2·`COARSE_GRAINING_LENGTH`·t)/(√2·`WIDTH`)], where `WIDTH` is the keyword of `&REGION`. m is 0 inside the explicit solvent, voids between the molecules included, and 1 far from it. For a flat uniform film of density `DENSITY` the mask equals the erf plane of `&REGION` placed `OFFSET` above the Gibbs dividing surface of the film, so a planar setup carries over with the same `WIDTH`. The mask multiplies S_diel and S_ion, is rebuilt whenever a listed atom has moved, and contributes an analytic force on the listed atoms. If `LOWER_EDGE` and `UPPER_EDGE` are also given, the planar mask multiplies on top. Only `METHOD SOLVENT_EXCLUDED`. ρ̄ is built in reciprocal space (one FFT per rebuild, one more for the forces), so the cost does not grow with `COARSE_GRAINING_LENGTH`.

List the water oxygens and the ions of the film with weight 1. Do not list the hydrogens or the electrode atoms. On Pt(111) under three water bilayers with one K⁺ (88 atoms; eight static frames and 100 MD steps) this list kept the implicit solvent out of the gap between the electrode and the film: max(ε − 1) there stayed at or below 1.0e-4. Listing the electrode atoms as well, with weight `DENSITY` divided by the atomic density of the metal, moves the boundary on the bare side of the slab outward (from 2.71 Å to 3.50 Å below the lowest Pt layer in that test) and lowered the work function of the charged slab by 0.29 eV. `ATOM_LIST` restricts the list to chosen atoms of the listed kinds.

| Keyword | Default | Meaning |
|---|---|---|
| `KIND <name> <weight>` | none | Repeatable. The film lists `Ow 1.0` and `K 1.0`. Kind name as in `&KIND` and its weight in ρ̄. At least one line is required; an unknown name aborts |
| `COARSE_GRAINING_LENGTH` | 2.4 Å | Gaussian width of the coarse graining; must exceed the `WIDTH` of `&REGION` |
| `DENSITY` | 0.0334 Å⁻³ | Number density of the explicit solvent in weighted atoms per volume (liquid water) |
| `OFFSET` | 2.4 Å | Distance of the boundary (m = ½) from the Gibbs dividing surface of the explicit solvent |
| `ATOM_LIST` | all atoms of the listed kinds | Restrict ρ̄ to these atoms (1-based); a list that selects no atom aborts |
| `FROM_DENSITY` | F | Build the field from the electron density instead: ρ̄ = g ∗ (1 − S_solv) with S_solv the solvent cavity function and `DENSITY` = 1. The mask then enters the Kohn–Sham potential, no atoms are listed, and `OFFSET` is counted from the surface S_solv = ½ with default 0 |

Each rebuild prints `SCCS| flexible region: atoms n  int rho_bar dV x  sum of weights y  int(1-m) dV [bohr^3] z`; x and y agree when the Gaussians are resolved on the grid, and z is the volume closed to the implicit solvent.

### `&SCCS &EXCLUSION`

Atom-following exclusion of the implicit dielectric and ions from explicit solvent atoms (explicit water, K⁺ or Na⁺ that move). Only `METHOD SOLVENT_EXCLUDED`; any other method aborts. M(r) = ∏ₐ[1 − ½ erfc((|r−Rₐ|−rₐ)/(√2·`WIDTH`))] multiplies S_diel and S_ion (the `&REGION` mask multiplies on top). The mask is rebuilt from the current positions at every SCCS call and contributes an analytic force on the excluded atoms. The cavitation surface is not masked. There are no default radii. Section parameter: `&EXCLUSION on`.

| Keyword | Default | Meaning |
|---|---|---|
| `KIND <name> <r_diel> <r_ion>` | none | Repeatable. The 203-atom Cu/water/K⁺ film uses `Ow 1.80 2.40`, `Hw 1.10 1.70`, `K 2.00 2.60`. Kind name as in `&KIND`, both radii in Å. At least one line is required; a kind without a line is not excluded; an unknown or repeated name aborts |
| `WIDTH` | 0.25 Å | Width of the erfc edge; must be > 0. The film uses 0.25 Å |
| `ATOM_LIST` | all atoms of the listed kinds | Restrict the exclusion to these atoms (1-based) |
| `PRINT_SUMMARY` | F | Per SCCS call: number of excluded atoms and ∫(1 − M_diel) dV, ∫(1 − M_ion) dV |

### Surface tension

`&SCCS GAMMA` keeps the CP2K default of 0. For solvation free energies of molecules set `GAMMA [mN/m] 18.47` with `ISODENSITY` or `GAMMA [mN/m] 11.09` with `SOLVENT_EXCLUDED` at `DIELECTRIC_RADIUS` 0.7 Å (1.1531 and 0.6923 meV/Å², 1 meV/Å² = 16.0218 mN/m), fitted on alkanes with the cavity settings of this document. The fit read the cavitation term from a run at `GAMMA [mN/m] 1.0` and scaled it; setting the value in the input also feeds the cavitation potential back into the SCF, a second-order difference not measured separately. For slabs, adsorption energies, potentials of zero charge and capacitances leave it at 0: at `GAMMA [mN/m] 14.50` the term moved the H adsorption energy on Pt(111) by only −1.8 meV.

### `&DFT &SCF &GCE`

CP2K's own keywords `TARGET_WORKFUNCTION` and `MIXING_COEF` belong to scheme `MU_MIXING`. The keywords below are the ones this code adds. The counter charge of the electrode is either `PLANAR_COUNTER_CHARGE` or an SCCS electrolyte (`DEBYE_LENGTH` > 0); both together abort. A constant-potential run needs `&DIAGONALIZATION`, `&SMEAR METHOD FERMI_DIRAC` and `ADDED_MOS`. Orbital transformation is not supported.

#### Scheme and target

| Keyword | Default | Meaning |
|---|---|---|
| `UPDATE_SCHEME` | `MU_MIXING` | `MU_MIXING`: Fermi level mixed toward the target work function. `CAPACITANCE`: N_e ← N_e + C (WF − WF_target) with a secant-updated capacitance C, Fermi–Dirac occupations at fixed N_e. `DYNAMICS`: electron-number dynamics, see below |
| `UPDATE_EVERY` | `SCF_STEP` | `CAPACITANCE`: update N_e at every SCF step, or once per converged SCF (`SCF_CYCLE`, for MD and `GEO_OPT`). `SCF_CYCLE` is the supported mode and the only one for `DYNAMICS` and `ENFORCE_TARGET` |
| `TARGET_POTENTIAL` | unset | Electrode potential in V vs SHE, entered as eV; sets the target work function to `PHI_SHE` + value. `CAPACITANCE` and `DYNAMICS` need one of `TARGET_POTENTIAL`, `TARGET_WORKFUNCTION` and `TARGET_CHEMICAL_POTENTIAL`, except for a fixed-N_e block (`MAX_OUTER_UPDATES` 0 with `UPDATE_EVERY SCF_CYCLE`) |
| `PHI_SHE` | 4.600 eV | Absolute SHE potential used with `TARGET_POTENTIAL`. 4.600 eV is the value calibrated for this implementation against the experimental pzc of nine Cu, Ag and Au faces (fit 4.602 eV) |
| `TARGET_CHEMICAL_POTENTIAL` | unset | Absolute electron chemical potential [eV], μ_e = −WF; overrides `TARGET_WORKFUNCTION` and `TARGET_POTENTIAL` |
| `EXCESS_ELECTRONS` | 0 | Fractional excess electrons added to the `CHARGE` state at the start; with `MAX_OUTER_UPDATES 0` this gives a constant-charge run at a fractional electron number |
| `SOLVENT_ONLY_FOR_CHARGE` | F | The implicit solvent only sets N_e; once N_e is fixed on an update step the solvent is switched off and the SCF re-converges in vacuum. On steps without an update (`UPDATE_INTERVAL`) it stays off, so energies and forces are solvent-free |

#### Update loop (`CAPACITANCE`, `UPDATE_EVERY SCF_CYCLE`)

| Keyword | Default | Meaning |
|---|---|---|
| `MAX_OUTER_UPDATES` | 20 | Cap on N_e updates per update step, each applied after a converged inner SCF. `MAX_SCF` counts the inner steps since the last update (and from the start of the SCF for the first N_e), so every update cycle gets the full `MAX_SCF` |
| `UPDATE_INTERVAL` | 1 (recommended for MD: 5, a balance between speed and accuracy) | Update only on force evaluations 1, n, 2n, 3n, …, counting the initial one as 1. For n > 1 the first update of an update evaluation starts from the last accepted capacitance and ignores the secant point of earlier ionic steps (it contains the geometry drift of the steps in between); the secant is rebuilt inside the evaluation. An evaluation that makes a single update (`MAX_OUTER_UPDATES` 1, `CONVERGENCE_MODE` 2, or mode 1 after the first step) uses the pair of the earlier update evaluation, with the new capacitance bounded to 0.5 to 2 times the last accepted one |
| `CONVERGENCE_MODE` | 3 | 1: converge the work function on the first ionic step only, afterwards one update per update step; 2: always one update per update step; 3: converge it on every update step (up to `MAX_OUTER_UPDATES`). `MAX_OUTER_UPDATES 0` keeps N_e fixed in every mode. Ignored by `DYNAMICS` |
| `FERMI_TOL` | 0.01 eV | Tolerance on \|WF − WF_target\| for the update loop, `CONVERGENCE_CHECK` and `ENFORCE_TARGET` |
| `CONVERGENCE_CHECK` | F | Also require \|WF − WF_target\| < `FERMI_TOL` for SCF convergence |
| `FINAL_TOL_FRAC` | 1 | After a converged inner SCF keep updating N_e while \|WF − WF_target\| > `FERMI_TOL` × value. Range (0, 1] |
| `CAPACITANCE_INIT` / `_MIN` / `_MAX` | 1 / 0.05 / 2 e/V | Initial value and clamp of the secant capacitance; all positive, MIN ≤ MAX |
| `NE_STEP_MAX` | 0.5 e | Largest N_e change per update |
| `FIRST_UPDATE_SIZE` | 0 = off | Fixed size [e] of the first update, sign from the work-function error; 0 takes the step from `CAPACITANCE_INIT` |
| `FERMI_UPDATE_TOL` | 0.02 eV | Refresh the secant capacitance only when the work function moved more than this since the last update |
| `CAP_REFRESH_MIN` | 0 eV | The secant is refreshed only when the work function moved more than max(`FERMI_UPDATE_TOL`, value) |
| `MIN_STEPS_AFTER_UPDATE` | 12 | Inner SCF steps that must follow an update before the SCF may be declared converged. With `PIN_MIXER_G0` off and Broyden mixing it is raised to `NBROYDEN` + 2 with a warning |
| `UPDATE_SCF_DELTA` | 0.01 | `UPDATE_EVERY SCF_STEP` only: update N_e only while the SCF density change of the previous step is below this |
| `PIN_MIXER_G0` | F; T for `DYNAMICS` unless given | Also set the G = 0 coefficient of the mixer input to N_e/V, so a new N_e takes effect in the density at once. Spin-restricted runs only; switched off with a warning for UKS |
| `RESPONSE_GUESS` | F | At an N_e update the density that enters the next inner SCF is ρ + dN_e·∂ρ/∂N_e, with ∂ρ/∂N_e the difference quotient of the densities at the last two updates of the same force evaluation (kept for later force evaluations). The mixer input is set to the same density. Until two such states exist the update only shifts the G = 0 coefficient. Spin-restricted runs only. It acts only when one force evaluation makes two or more N_e updates: in a 323-atom Cu(111)/water film it lowered the SCF step count of the first force evaluation from 163 to 138, and to 74 together with early updates (122 with early updates alone). Smaller test cells that make one update per force evaluation were unchanged, and early updates forced on them cost SCF steps. Each use prints `GCE\| response guess: density moved along d rho/d N_e by x e` |
| `KEEP_MIXER_BELOW` | 0; 0.05 e for `DYNAMICS` unless given | An N_e update with \|dN_e\| [e] below this keeps the Broyden history instead of restarting it |
| `EARLY_UPDATE_DELTA` | 0 = off | If > 0, N_e is also updated before the inner SCF has converged, once the density change is below this and `EARLY_UPDATE_MIN_STEPS` steps have run at the current N_e. The final N_e is still judged on a converged density |
| `EARLY_UPDATE_MIN_STEPS` | 3 | Steps at the new N_e before an early update is allowed |
| `EARLY_UPDATE_MAX` | 2 | Largest number of early updates per SCF run; they do not count against `MAX_OUTER_UPDATES` |
| `EARLY_UPDATE_WF_STABLE` | 0 eV | If > 0, an early update also needs the work function to change by less than this between two inner steps |
| `EARLY_UPDATE_SECANT` | F | The first early update records its (N_e, WF) pair so that the next update has a secant capacitance, and later early updates refresh it. `F`: early updates leave the secant history to the first converged update |
| `CARRY_GAIN` | 0 = off | After each converged ionic step N_e moves by value × C × (WF − WF_target) for the next step. 0.5 to 1 are sensible. Not combinable with `DYNAMICS` |
| `ENFORCE_TARGET` | F | An update step that used all `MAX_OUTER_UPDATES` and is still off the target by more than `FERMI_TOL` stops the run with the reason (target inside a gap, steps too small). In `GEO_OPT` every force evaluation is a full update step, so the last one is within `FERMI_TOL`. No effect with `UPDATE_EVERY SCF_STEP` or `MAX_OUTER_UPDATES` 0 |

#### Restart state

The code writes these after every ionic step into the `&GCE` section of the running input, so every MD and `GEO_OPT` restart file carries them. They are not meant to be edited.

| Keyword | Default | Content |
|---|---|---|
| `RESTART_VALID` | F | The `RESTART_*` values are valid (shared by `CAPACITANCE` and `DYNAMICS`) |
| `RESTART_N_E`, `RESTART_N_INPUT`, `RESTART_CAPACITANCE`, `RESTART_ION_STEP` | 0 | N_e [e] of the last evaluation, N_e of the input charge state [e], capacitance, number of completed force evaluations (shared by both schemes) |
| `RESTART_N_LAST`, `RESTART_WF_LAST`, `RESTART_SECANT_READY`, `RESTART_LAST_WF`, `RESTART_STEPS_SINCE_UPDATE` | 0 / F | `CAPACITANCE`: previous secant point (N_e, WF), secant flag, work function of the last ionic step, inner steps since the last update |
| `RESTART_P_N`, `RESTART_DU_LAST`, `RESTART_NEXT_LEN`, `RESTART_OPEN_PENDING`, `RESTART_XI`, `RESTART_P_XI` | 0 / 1 / F | `DYNAMICS`: momentum of N_e, last U − U_target, MD steps to the next potential measurement, opening half step due, chain positions and momenta |

`&EXT_RESTART` has a keyword `RESTART_GCE` (like the other `RESTART_*` keywords of that section it follows `RESTART_DEFAULT`, which is T) that copies the explicit `RESTART_*` keywords of the restart file into the new input.

### `&DFT &SCF &GCE UPDATE_SCHEME DYNAMICS`

`RUN_TYPE MD` only, `UPDATE_EVERY SCF_CYCLE`, fixed MD time step. The electron number n is a dynamical variable (mass M_n, momentum P_n) driven by U − U_T, one SCF per ionic step and no N_e update inside the SCF after the first force evaluation. The first evaluation converges n to the target with the usual secant loop (`MAX_OUTER_UPDATES`; 0 starts from `EXCESS_ELECTRONS`):

```
dn/dt = P_n/M_n      dP_n/dt = (U − U_T) − P_n P_ξ1/Q_1      (chain: dξ_j/dt = P_ξj/Q_j, dP_ξj/dt = P_ξ(j−1)²/Q_(j−1) − kT − P_ξj P_ξ(j+1)/Q_(j+1))
M_n = 1/(C ω²),  ω = 2π/POTENTIOSTAT_PERIOD,  Q_j = kT/ω_th²
```

The chain uses the symmetric Suzuki–Yoshida decomposition of https://doi.org/10.1080/00268979600100761. The MD time step and temperature are read from `MOTION&MD`.

| Keyword | Default | Meaning |
|---|---|---|
| `POTENTIOSTAT_PERIOD` | required, > 0 fs | Period 2π/ω of the free oscillation of n; sets M_n. Choose it well above 10 MD time steps (the code warns below 10 update cycles) |
| `POTENTIOSTAT_CAPACITANCE` | 0 = secant value measured by the first evaluation | C = dN_e/dU of the cell [e/V]; the run stops if none is available. Measure it with a static N_e scan for production |
| `POTENTIOSTAT_CHAIN_LENGTH` | 0 | 0: constant energy (K_n + Ω conserved); 1–16: Nosé–Hoover chain at `POTENTIOSTAT_TEMPERATURE` (3 or more for statistics; a single link does not sample a harmonic mode ergodically) |
| `POTENTIOSTAT_TEMPERATURE` | `MOTION&MD TEMPERATURE` | kT of the chain |
| `POTENTIOSTAT_THERMOSTAT_PERIOD` | `POTENTIOSTAT_PERIOD` | Q_j = kT/ω_th² |
| `POTENTIOSTAT_INITIAL_RATE` | 0 e/fs | Initial dn/dt; a continuation run restores the state of the previous run instead |

With `SOLVENT_ONLY_FOR_CHARGE` the potential is measured on the steps with the solvent on; on those steps the `GCE| trace` line, `RESTART_LAST_WF` and `CARRY_GAIN` use the solvent-on work function (the steps without solvent have no such value and report the solvent-off one, which is not the electrode potential); `UPDATE_INTERVAL` n sets the spacing and kicks are applied as impulses of n steps. The conserved quantity is (MD `Cons Qty`) + `Hextra`, with Ω = E + WF_T (n − n_input).

### Restart of constant-potential runs

Continue with the `.restart` file written by the MD driver as the input (`cp2k -i <project>-1.restart`), or from a fixed deck with `&EXT_RESTART`. The first force evaluation repeats the last evaluation of the previous run at the same N_e and the same update-step phase, no extra N_e update is made, and the next steps follow the uninterrupted run. If the last evaluation of the previous run ended off the target (all `MAX_OUTER_UPDATES` used), the repeated evaluation makes one more update. The restored N_e is only a starting value: the update steps and `CARRY_GAIN` compare the work function with the `TARGET_POTENTIAL` of the new input, so after a restart with another target potential N_e moves to it within `UPDATE_INTERVAL` force evaluations. `RESTART_GCE F` skips the electron-number state and N_e starts again from `CHARGE` and `EXCESS_ELECTRONS`. A restart for another `CHARGE` or atom count aborts with a message. With `DYNAMICS`, `POTENTIOSTAT_CHAIN_LENGTH` must equal the length stored in the restart file.

### Spin

With `UPDATE_SCHEME CAPACITANCE` or `DYNAMICS`, `UKS` shares one electron number and one Fermi level between the spins. `FIXED_MAGNETIC_MOMENT` with a constant-potential run aborts. With the common Fermi level `MULTIPLICITY` only sets the initial occupation (a warning says so): a run started from a triplet can converge to the non-magnetic state, and a fixed-N_e run from the same triplet can reach a high-spin solution that lies lower in energy. Scheme `MU_MIXING` with `LSD` keeps the per-spin occupation call.

### `&DFT &PRINT &SCCS`

`IONIC_CHARGE_DENSITY` (debug print level) writes cube files of the ionic charge density of the Poisson–Boltzmann electrolyte (`DEBYE_LENGTH` > 0), with `STRIDE` and `APPEND`.


## Former keyword names

The keywords below were renamed. The former names still work in every input and restart file; the code reads either name and writes the new one.

| Section | Former name | Name |
|---|---|---|
| `&SCCS` | `ION_CHARGE` | `ION_VALENCE` |
| `&SCCS` | `EPS_NLPB` | `SOLVER_TOLERANCE` |
| `&SCCS` | `MAX_ITER_NLPB` | `SOLVER_MAX_ITER` |
| `&SCCS` | `NLPB_ENERGY_TOL` | `SOLVER_ENERGY_TOL` |
| `&SCCS` | `MAX_NLPB_FAILURES` | `SOLVER_MAX_FAILURES` |
| `&SCCS` | `DEBUG_NLPB` | `SOLVER_DEBUG` |
| `&SOLVENT_EXCLUDED` | `RADIUS_SOLVENT` | `SOLVENT_RADIUS` |
| `&SOLVENT_EXCLUDED` | `RADIUS_ION` | `ION_RADIUS` |
| `&SOLVENT_EXCLUDED` | `RADIUS_DIELECTRIC` | `DIELECTRIC_RADIUS` |
| `&SOLVENT_EXCLUDED` | `RADIUS_CAVITATION` | `CAVITATION_RADIUS` |
| `&REGION` | `Z_LOW` | `LOWER_EDGE` |
| `&REGION` | `Z_HIGH` | `UPPER_EDGE` |
| `&REGION &FLEXIBLE` | `SIGMA` | `COARSE_GRAINING_LENGTH` |
| `&GCE` | `NELECTRON_EXCESS` | `EXCESS_ELECTRONS` |

`ION_RADIUS` written directly under `&SCCS` remains an alias of `ION_STERIC_RADIUS` (the lattice radius of the size-modified ions). Inside `&SOLVENT_EXCLUDED` it is the closest approach of an ion to the solute surface.

## Symbols

Symbol, quantity, input keyword and unit. Quantities without an input keyword are not listed.

| Symbol | Quantity | Keyword | Unit |
|---|---|---|---|
| **`&SCCS`** | | | |
| ε_b | Bulk relative permittivity | `RELATIVE_PERMITTIVITY` | – |
| ε_∞ | Optical permittivity of the solvent | `EPS_INF` | – |
| λ_D | Debye length | `DEBYE_LENGTH` | Å |
| c | Ion concentration | `ION_CONCENTRATION` | M |
| z_ion | Ion valence (symmetric z_ion:z_ion electrolyte) | `ION_VALENCE` | – |
| T | Electrolyte temperature | `TEMPERATURE` | K |
| R_h | Hydrated-ion radius of the lattice gas | `ION_STERIC_RADIUS` | Å |
| ψ | Symmetrised form of the Poisson operator | `PSI_FORM` | – |
| τ | Cavitation surface tension | `GAMMA` | mN/m |
| **`&SCCS ISODENSITY and SOLVENT_EXCLUDED`** | | | |
| n_c | Electron-density isovalue of the solute surface | `RHO_CAVITY` | Å⁻³ |
| σ_c | Width of the solute surface on the ln(n/n_c) axis | `SWITCH_WIDTH` | – |
| **`&SCCS SOLVENT_EXCLUDED`** | | | |
| λ_p | Decay length of the probe-shell kernel | `PROBE_KERNEL_WIDTH` | Å |
| R_solv | Closest approach of a solvent molecule centre to the solute surface | `SOLVENT_RADIUS` | Å |
| R_ion | Closest approach of an electrolyte ion to the solute surface | `ION_RADIUS` | Å |
| R_diel | Inward extension of dielectric screening from the solvent-centre surface | `DIELECTRIC_RADIUS` | Å |
| R_cav | Inward displacement of the surface that carries the cavitation term | `CAVITATION_RADIUS` | Å |
| **`&SCCS EXCLUSION`** | | | |
| w_ex | Edge width of the exclusion sphere | `WIDTH` | Å |
| **`&SCCS REGION`** | | | |
| z_low, z_high | Lower and upper edge of the implicit region | `LOWER_EDGE`, `UPPER_EDGE` | Å |
| w | Edge width of the region | `WIDTH` | Å |
| **`&SCCS REGION FLEXIBLE`** | | | |
| ℓ | Gaussian width of the coarse-grained solvent density | `COARSE_GRAINING_LENGTH` | Å |
| ρ_0 | Number density of the explicit solvent | `DENSITY` | Å⁻³ |
| d_off | Distance of the boundary from the Gibbs dividing surface | `OFFSET` | Å |
| **`&DFT SCF GCE`** | | | |
| Φ_SHE | Absolute potential of the SHE | `PHI_SHE` | eV |
| U_T | Target electrode potential vs SHE | `TARGET_POTENTIAL` | eV (as V) |
| μ_e | Target electron chemical potential | `TARGET_CHEMICAL_POTENTIAL` | eV |
| N_exc | Excess electrons relative to the neutral cell | `EXCESS_ELECTRONS` | e |
| δ_W | Tolerance on |W − W_T| | `FERMI_TOL` | eV |
| C_0 | Initial differential capacitance | `CAPACITANCE_INIT` | e/V |
| C_min, C_max | Lower and upper clamp of the secant capacitance | `CAPACITANCE_MIN`, `CAPACITANCE_MAX` | e/V |
| ΔN_max | Largest change of N_e per update | `NE_STEP_MAX` | e |
| k_upd | Ionic steps between outer updates | `UPDATE_INTERVAL` | steps |
| n_out | Maximum number of outer N_e updates per active step | `MAX_OUTER_UPDATES` | – |
| g_ff | Feed-forward gain of the electron-number correction carried to the next ionic step | `CARRY_GAIN` | – |
| τ_pot | Potentiostat period | `POTENTIOSTAT_PERIOD` | fs |
| **`PLANAR_COUNTER_CHARGE (CP2K)`** | | | |
| d_pl | Distance of a counter-charge plane from the cell boundary | `DIST_EDGE` | Å |
| σ_pl | Gaussian width of the planar counter charge | `GAU_C` | Å |

