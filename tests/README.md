# Regression tests

Seven directories in the CP2K layout, read through `TEST_FILES.toml`.

| Directory | What it pins | Cost |
|---|---|---|
| `regtest-sccs-lpb-1` | Energies on H2O (10 Å box, 400 Ry) for both cavity models, both solvers, the linear electrolyte, size-modified ions with the nonlinear dielectric, the pseudo-core keywords and the `&DEVEL` section | seconds to minutes |
| `regtest-sccs-lpb-2` | `RUN_TYPE DEBUG` finite-difference forces (600 Ry, 8 Å box): the pseudo-core force term, the linear electrolyte, and the chained cavity derivatives of `SOLVENT_EXCLUDED` | minutes; above the usual dashboard budget |
| `regtest-sccs-lpb-3` | Constant potential on a 12-atom Al(111) slab: the outer loop runs, prints `GCE| outer update`, and the energy is reproducible | minutes |
| `regtest-sccs-lpb-4` | Guards: a neutral GCE block on a k-point mesh runs, `SOLVER_MAX_FAILURES` aborts after consecutive NLPB failures (`H2O_nlpb_failure_abort.inp` must abort, run by hand) | minutes |
| `regtest-sccs-lpb-5` | `&SCCS&EXCLUSION`: K+ + 6 H2O `RUN_TYPE DEBUG` forces (linear PB; nonlinear dielectric + nonlinear ions; fully linear branch), H2O `RUN_TYPE DEBUG` with large spheres (pins the explicit mask force), EXCLUSION OFF energy identical to the base deck, EXCLUSION with ISODENSITY aborts | DEBUG decks: tens of minutes (19 atoms x 6 displacements, CUTOFF 600) |
| `regtest-sccs-lpb-6` | GCE state and dynamics: restart state (MD and `EXT_RESTART`), the per-step `GCE| omega` and `GCE| conserved` lines, `ENFORCE_TARGET` (unreachable target, geometry optimisation), UKS | MD decks: 12 steps of a 12-atom slab, 3 min on 16 ranks (NERSC, Ibex 5 min); the others 1–3 min on 8 ranks |
| `regtest-sccs-lpb-7` | `&SCCS&REGION&FLEXIBLE` on H2O (600 Ry, 8 Å box): energy with the mask built from the O atom, `RUN_TYPE DEBUG` forces including the mask force, three input guards that abort | energy 5 min on 2 ranks; DEBUG 19 energy evaluations |

## Former keyword names and default PHI_SHE

Three decks end in `_legacy_names`: `regtest-sccs-lpb-1/H2O_solvent_excluded_nli_nld_region_probe_legacy_names.inp`, `regtest-sccs-lpb-6/Al111_gce_fixq_mode2_legacy_names.inp` and `regtest-sccs-lpb-7/H2O_solvent_excluded_flexible_legacy_names.inp`. They are the decks of the same name without the suffix with every renamed keyword written under its former name (together all 14), and give the same energy. `regtest-sccs-lpb-6/Al111_gce_rks_default_phi.inp` is `Al111_gce_rks.inp` without `PHI_SHE`: the target work function uses the default 4.600 eV.

## Running

Copy the seven directories to `tests/QS/` of the CP2K tree, add them to `tests/TEST_DIRS`, and run the CP2K regression tester as usual.

`regtest-sccs-lpb-3` has a second condition that a matcher cannot express: the output must contain at least one `GCE| outer update` line.

```bash
grep -c "GCE| outer update" Al111_gce_cap_lpb.inp.out    # must be > 0
```

`regtest-sccs-lpb-4` holds one deck whose expected result is an abort. Run it by hand and check the exit code and the message (header of the deck):

```bash
cp2k.psmp -i H2O_nlpb_failure_abort.inp -o H2O_nlpb_failure_abort.out; echo $?     # non-zero
grep -c "SOLVER_MAX_FAILURES" H2O_nlpb_failure_abort.out                             # > 0
```

`regtest-sccs-lpb-7` holds three decks whose expected result is an abort. Run them by hand; the message is in the header of each deck:

```bash
cp2k.psmp -i H2O_flexible_sigma_abort.inp -o a.out; echo $?; grep -c "COARSE_GRAINING_LENGTH must exceed the REGION WIDTH" a.out   # non-zero, > 0
cp2k.psmp -i H2O_flexible_kind_abort.inp -o b.out; echo $?; grep -c "unknown atomic kind name XX" b.out           # non-zero, > 0
cp2k.psmp -i H2O_flexible_zlow_abort.inp -o c.out; echo $?; grep -c "needs both LOWER_EDGE and UPPER_EDGE" c.out           # non-zero, > 0
```

`regtest-sccs-lpb-6` (GCE state and dynamics): the energies are in `TEST_FILES.toml`; the conditions that a matcher cannot express and the decks that are run by hand (two aborts, three restarts from `Al111_gce_md_state-1_7.restart`) are listed in its header. The restart decks need the restart file of `Al111_gce_md_state.inp` in their directory; the `RESTART_*` lines of that file and the `GCE| state read from the restart` line of the restarted run hold the same values, and the electron numbers of the steps after the repeated one follow the uninterrupted run.

```bash
grep -c "GCE| omega" Al111_gce_md_state.out; grep -c "GCE| conserved" Al111_gce_md_state.out   # 13 and 13
cp2k.psmp -i Al111_gce_unreach_abort.inp -o abort.out; echo $?                                 # non-zero
grep -c "stays on one side" abort.out                                                           # > 0
cp2k.psmp -i H2O_gce_gap_abort.inp -o gap.out; echo $?                                         # non-zero
grep -c "crosses the target without hitting it" gap.out                                         # > 0
```

`regtest-sccs-lpb-1/H2O_solvent_excluded_nli_nld_region_probe.inp` (PROBE_ATOMS, PRINT_POL_DIAGNOSTICS) reproduces the `region_apply` energy; `grep -c "SCCS| region post-mask"` on its output must be > 0. `Al111_gce_kpoints_neutral.inp` has no reference energy; it is registered with an empty matcher list, so the tester checks that it runs to completion.

`regtest-sccs-lpb-5`: `H2O_se_exclusion_off.inp` must give the energy of the base deck (-17.23505472278111). The three `K6H2O_dbg_*` decks compare analytic and finite-difference forces on all 19 atoms, which are all excluded atoms; criterion per atom <= 1E-4 Ha/bohr and <= 0.1 % relative (the M042 sum is only the dashboard handle). Run each DEBUG deck on 1 and on 32 MPI ranks: the forces must agree to 1E-8 Ha/bohr (this checks that the exclusion force is not double counted over ranks). `H2O_isodensity_exclusion_abort.inp` must abort with `SCCS EXCLUSION is implemented only for METHOD SOLVENT_EXCLUDED`:

```bash
cp2k.psmp -i H2O_isodensity_exclusion_abort.inp -o abort.out; echo $?   # non-zero
grep -c "EXCLUSION is implemented only" abort.out                       # > 0
```

Reference values of the DEBUG decks: to be generated on HPC. The radii in these decks (Ow 1.8/2.4, Hw 1.1/1.7, K 2.0/2.6 A) are test values, not calibrated.

## Decks that are not listed in `TEST_FILES.toml`

The CP2K tester counts a non-zero exit status as a failure, so decks that must abort cannot be listed. All other decks of the six directories are listed except the following.

| Deck | Why it is not listed |
|---|---|
| `regtest-sccs-lpb-4/H2O_nlpb_failure_abort.inp` | must abort |
| `regtest-sccs-lpb-5/H2O_isodensity_exclusion_abort.inp` | must abort |
| `regtest-sccs-lpb-6/Al111_gce_unreach_abort.inp`, `H2O_gce_gap_abort.inp`, `Al111_gce_dyn_{chain,energy,noperiod,scfstep}_abort.inp` | must abort |
| `regtest-sccs-lpb-6/Al111_gce_md_state_restart.inp`, `_restart_cold.inp`, `_restart_newU.inp` | need the restart file of `Al111_gce_md_state.inp` and are checked by `grep` on the output |
| `regtest-sccs-lpb-5/K6H2O_dbg_se_exclusion_{lpb,nli_nld,linear_all}.inp` | `RUN_TYPE DEBUG` decks of tens of minutes without a recorded `M042` reference |

## Tolerances

- Linear and fixed-point energy paths: relative 1.0E-10.
- Newton paths with a nonlinear response: relative 1.0E-7.
- Finite-difference force residual (`M042`, the sum over 3N components of |F_numerical − F_analytical|): 1 %. It pins a physical residual of order 1E-7 Ha/bohr; a missing force term moves it by orders of magnitude.
- Constant-potential energy: 1.0E-4, because a metallic grand-canonical SCF scatters by about 1E-4 Ha between builds.

## Known weaknesses of the suite

- `H2O_isodensity_core_density.inp` (e) is numerically identical to `H2O_isodensity_jacobi.inp` (a): the two keywords it sets explicitly are now the defaults. It checks parsing only. Isolating the effect of the pseudo-core on the cavity needs a companion deck with `ADD_CORE_DENSITY F`.
- `H2O_isodensity_jacobi.inp` (a) and `H2O_isodensity_newton.inp` (b) differ by 2.0E-6 Ha because the Jacobi solver stops on the step-to-step change. Each pins its own value.
- There is no `RUN_TYPE DEBUG` force test for the constant-potential path.

## Five-point charging validation (`tests/fivepoint/`)

The one physics test that must reproduce after any renaming or refactoring: Pt(111) 4×4×4 (64 atoms, 106.4615 Å² per cell),
`METHOD SOLVENT_EXCLUDED`, 1 M size-modified ions (`ION_STERIC_RADIUS 4.0 Å`), nonlinear dielectric, `ADD_CORE_DENSITY F`,
`&GCE` with `EXCESS_ELECTRONS` −1, −0.5, 0, +0.5, +1 e. `fit_fivepoint.py fp_*.out` converts the Fermi energies to
U = −E_F − 4.43 V vs SHE, fits N_exc = a·U + b and reports C = −a, U_pzc = b/C and ΔΩ(U) = −½C(U − U_pzc)².

Reference (24 MPI ranks × 2 OpenMP per point, 6–8 min each):

| N_exc [e] | E_F [Ha] |
|---|---|
| −1.0 | −0.22556788447 |
| −0.5 | −0.20937012078 |
| 0.0 | −0.19375886315 |
| +0.5 | −0.17863256591 |
| +1.0 | −0.16284938159 |

C = 1.1764 e/V = 17.70 µF/cm², U_pzc = 0.8500 V vs SHE, R² = 0.999880. Separate MPI runs reproduce E_F to ~1e-11 Ha;
the script passes at |ΔE_F| < 1e-7 Ha. Run each deck with `cp2k.psmp -i fp_x00.inp -o fp_x00.out` in this directory
(`cell.inc` and `pt111_4x4x4.coord` are included).
